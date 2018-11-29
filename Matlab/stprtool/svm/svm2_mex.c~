/*---------------------------------------------------------------------------
 svm2_mex.c: MEX-file for binary SVM with L2-soft margin solver.

 Compile: 
  mex svm2_mex.c gnppsolver.c kernel_fun.c

 Synopsis:
  [Alpha,bias,exitflag,kercnt,access,trnerr,t,UB,LB,History] = 
     svm2_mex(data,labels,ker,arg,C,solver,tmax,tolabs,tolrel,thlb,cache,verb)

 Input:
  data [dim x num_data] Training vectors.
  labels [1 x num_data] Labels.
  ker [string] Kernel identifier.
  arg [1 x nargs] Kernel argument.
  C [1x1] Regularization constant.
  solver [string] Solver; options are 'mdm'.
  tmax [1x1] Maximal number of iterations.
  tolabs [1x1] Absolute tolerance stopping condition.
  tolrel [1x1] Relaitve tolerance stopping condition.
  thlb [1x1] Threshold on lower bound.
  cache [1x1] Number of columns of kernel matrix to be cached.
    It takes cache*num_data*size(double) bytes of memory.
  verb [1x1] If 1 then some info about the training is printed.

 Output:
  Alpha [nclass x num_data] Weights.
  bias [1x1] Bias.
  exitflag [1x1] Indicates which stopping condition was used:
    UB-LB <= tolabs           ->  exit_flag = 1   Abs. tolerance.
    (UB-LB)/(LB+1) <= tolrel  ->  exit_flag = 2   Relative tolerance.
    t >= tmax                 ->  exit_flag = 0   Number of iterations.
  kercnt [1x1] Number of kernel evaluations.
  access [1x1] Number or requested columns of the kernel matrix.
  trnerr [1x1] Training error.
  t [1x1] Number of iterations.
  UB [1x1] Upper bound on the optimal solution.
  LB [1x1] Lower bound on the optimal solution.
  History [2x(t+1)] UB and LB with respect to number of iterations.

 Modifications:
 27-nov-2004, VF
-------------------------------------------------------------------- */

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "kernel_fun.h"

#define INDEX(ROW,COL,DIM) ((COL*DIM)+ROW)

#define MINUS_INF INT_MIN
#define PLUS_INF  INT_MAX

#define KDELTA(A,B) (A==B)
#define KDELTA4(A1,A2,A3,A4) ((A1==A2)||(A1==A3)||(A1==A4)||(A2==A3)||(A2==A4)||(A3==A4))


/* Declaration of global variables */

unsigned long access_cnt;  /* counter of access to H matrix */
long num_data;       /* number of input training examples */
double reg_const;    /* regularization constant */
double *vector_y;    /* Pointer to labels */

long Cache_Size;     /* number of cached columns (min 1) */

/* cache (FIFO) for columns of the kernel matrix */
long *cache_index;                  /* indices cached of kernel columns */
long first_kernel_inx;              /* index of first inserted column */
double **kernel_columns;            /* pointers at cached columns */

/* ------------------------------------------------------------
  Returns pointer at a-th column of the matrix H = (y*y')*K.
  This function maintains FIFO cache of kernel columns.

  (note: the b-th column must be preserved in the cache during 
   updating but b is from (a(t-2), a(t-1)) where a=a(t) and
   thus FIFO with more than three columns does not have to 
   take care od b.)
------------------------------------------------------------ */
void *get_col( long a, long b ) 
{
  double *col_ptr;
  double y;
  long i;
  long inx;

  access_cnt = access_cnt + 1;

  inx = -1;
  for( i=0; i < Cache_Size; i++ ) {
    if( cache_index[i] == a ) { inx = i; break; }
  }
    
  if( inx != -1 ) {
    col_ptr = kernel_columns[inx];
    return( col_ptr );
  }
   
  col_ptr = kernel_columns[first_kernel_inx];
  cache_index[first_kernel_inx] = a;

  first_kernel_inx++;
  if( first_kernel_inx >= Cache_Size ) first_kernel_inx = 0;

  y = vector_y[a];
  for( i=0; i < num_data; i++ ) {
    if( vector_y[i] == y )  
    {
      col_ptr[i] = 2*kernel(i,a); 
    }
    else 
    {
      col_ptr[i] = -2*kernel(i,a);
    }
  }

  col_ptr[a] = col_ptr[a] + reg_const;

  return( col_ptr );
}

/* -------------------------------------------------------------------
 Main MEX function - interface to Matlab.
-------------------------------------------------------------------- */
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[] )
{
  char solver[20];   /* solver identifier */
  int exitflag;      /* output arg */
  int buf_len;       /* real length of the solver identifier */
  long i;            /* common use loop variables */
  long tmax;         /* input arg - max number of iteration */ 
  long t;            /* output arg - number of iterations */
  long verb;         /* input argument */
  double nconst;        
  double b;          /* bias */
  double C;          /* input arg - regularization const */
  double thlb;       /* input arg */
  double tolrel;     /* input arg */
  double tolabs;     /* input arg */
  double aHa11;
  double aHa22;
  double trnerr;     /* output arg */
  double *vector_c;  /* auxiliary */ 
  double *alpha;     /* solution vector */ 
  double *History;   /* output arg */
  double *diagK;     /* cache for diagonal of virtual K matrix */
  double *tmp_ptr;

  /*------------------------------------------------------------------- */
  /* Take input arguments                                               */
  /*------------------------------------------------------------------- */

  if( nrhs != 12) mexErrMsgTxt("Incorrect number of input arguments.");

  dataA = mxGetPr(prhs[0]);   /* pointers at data */
  dataB = dataA;
  dim = mxGetM(prhs[0]);      /* data dimension */
  num_data = mxGetN(prhs[0]); /* number of data */
  vector_y = mxGetPr(prhs[1]);  /* pointer at data labels */

  /* take kernel identifier and its argument */
  ker = kernel_id( prhs[2] ); 
  if( ker == -1 ) mexErrMsgTxt("Improper kernel identifier.");
  arg1 = mxGetPr(prhs[3]);

  C = mxGetScalar(prhs[4]);   /* regularization constant */

  /* take string identifier QP solver to be used */
  if( mxIsChar( prhs[5] ) != 1) mexErrMsgTxt("solver must be string.");
  buf_len = (mxGetM(prhs[5]) * mxGetN(prhs[5])) + 1;
  buf_len = (buf_len > 20) ? 20 : buf_len;
  mxGetString( prhs[5], solver, buf_len );

  /* maximal allowed number of iterations */
  tmax = mxIsInf( mxGetScalar(prhs[6])) ? INT_MAX : (long)mxGetScalar(prhs[6]); 
  tolabs = mxGetScalar(prhs[7]);   /* abs. precision defining stopping cond*/
  tolrel = mxGetScalar(prhs[8]);   /* rel. precision defining stopping cond*/
  /* threshold on lower bound */
  thlb = mxIsInf( mxGetScalar(prhs[9])) ? DBL_MAX : (double)mxGetScalar(prhs[9]); 

  Cache_Size = (long)mxGetScalar(prhs[10]);  /* cache size */
  if( Cache_Size < 3 ) mexErrMsgTxt("Cache must be greater than 3."); 
  if( Cache_Size > num_data ) Cache_Size = num_data; 

  /* threshold on lower bound */
  verb = mxIsInf( mxGetScalar(prhs[11])) ? INT_MAX : (double)mxGetScalar(prhs[11]); 

  /*------------------------------------------------------------------- */
  /* Inicialization (caches, etc.)                                      */
  /*------------------------------------------------------------------- */

  /* constant added to diagonal of separable problem */
  if( C!=0 ) reg_const = 1/C; else reg_const = 0;

  ker_cnt = 0;            /* counter of kernel evaluations */
  access_cnt = 0;  /* counter for access to columns of kernel matrix */

  /* allocattes and precomputes diagonal of virtual K matrix */
  diagK = mxCalloc(num_data, sizeof(double));
  if( diagK == NULL ) mexErrMsgTxt("Not enough memory.");
  for(i = 0; i < num_data; i++ ) {
    diagK[i] = 2*kernel(i,i) + reg_const;
  }

  /* allocates memory for kernel cache */
  kernel_columns = mxCalloc(Cache_Size, sizeof(double*));
  if( kernel_columns == NULL ) mexErrMsgTxt("Not enough memory.");
  cache_index = mxCalloc(Cache_Size, sizeof(double));
  if( cache_index == NULL ) mexErrMsgTxt("Not enough memory.");

  for(i = 0; i < Cache_Size; i++ ) 
  {
    kernel_columns[i] = mxCalloc(num_data, sizeof(double));
    if(kernel_columns[i] == NULL) mexErrMsgTxt("Not enough memory.");

    cache_index[i] = -2;
  }
  first_kernel_inx = 0;

  /* Solution vector */
  plhs[0] = mxCreateDoubleMatrix(num_data,1,mxREAL);
  alpha = mxGetPr(plhs[0]);

  /* Vector c; for this problem set to zero */
  vector_c = mxCalloc(num_data, sizeof(double));
  if( vector_c == NULL ) mexErrMsgTxt("Not enough memory.");
  for(i = 0; i < num_data; i++ ) vector_c[i] = 0;


  /*------------------------------------------------------------------- */
  /* Call QP solver                                                     */
  /*------------------------------------------------------------------- */

  if ( strcmp( solver, "mdm" )==0 ) {  
     exitflag = gnpp_mdm( &get_col, diagK, vector_c, vector_y, num_data, 
         tmax, tolabs, tolrel, thlb, alpha, &t, &aHa11, &aHa22, 
         &History, verb ); 
  } else if ( strcmp( solver, "imdm" )==0 ) {  
     exitflag = gnpp_imdm( &get_col, diagK, vector_c, vector_y, num_data, 
         tmax, tolabs, tolrel, thlb, alpha, &t, &aHa11, &aHa22, 
         &History, verb ); 
  } else {
     mexErrMsgTxt("Unknown solver identifier.");
  }

  /*------------------------------------------------------------------- */
  /* Generate outputs                                                   */
  /*------------------------------------------------------------------- */

  /* computes alpha and b to be paramaters of the separating
     hyperplane in the cannonical form. */
  nconst = History[INDEX(1,t,2)];
  trnerr = 0; /* counter of training error */
  for(i = 0; i < num_data; i++ )
  {
     if(vector_y[i] == 1) 
     {
       alpha[i] = alpha[i]*2/nconst;
       if( alpha[i]/(2*C) >= 1 ) trnerr++;
     }
     else
     {
       alpha[i] = -alpha[i]*2/nconst;
       if( alpha[i]/(2*C) <= -1 ) trnerr++;
     }
  }
  
  /* bias vector b [1 x 1] */
  b = 0.5*(aHa22 - aHa11)/nconst;;
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[1])) = (double)b;

  /* exit_flag [1x1] */
  plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[2])) = (double)exitflag;

  /* kercnt [1x1] */
  plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[3])) = (double)ker_cnt;

  /* access [1x1] */
  plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[4])) = (double)access_cnt;

  /* trnerr [1x1] */
  plhs[5] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[5])) = trnerr;

  /* t [1x1] */
  plhs[6] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[6])) = (double)t;

  /* UB [1x1] */
  plhs[7] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[7])) = History[INDEX(1,t,2)];

  /* LB [1x1] */
  plhs[8] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[8])) = History[INDEX(0,t,2)];

  /* History [2 x (t+1)] */
  plhs[9] = mxCreateDoubleMatrix(2,t+1,mxREAL);
  tmp_ptr = mxGetPr( plhs[9] );
  for( i = 0; i <= t; i++ ) {
     tmp_ptr[INDEX(0,i,2)] = History[INDEX(0,i,2)];
     tmp_ptr[INDEX(1,i,2)] = History[INDEX(1,i,2)];
  }

  /*------------------------------------------------------------------- */
  /* Free used memory                                                   */
  /*------------------------------------------------------------------- */
  mxFree( vector_c );
  mxFree( History );
  mxFree( diagK );
  for(i = 0; i < Cache_Size; i++ ) mxFree(kernel_columns[i]);
  mxFree( kernel_columns );
  mxFree( cache_index );
}

