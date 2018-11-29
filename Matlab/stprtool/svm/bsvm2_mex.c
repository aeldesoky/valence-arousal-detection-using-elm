/*---------------------------------------------------------------------------
 bsvm2_mex.c: MEX-file for multi-class B-SVM with L2-soft margin.

 Compile: 
  mex new_bsvm2_mex.c gmnpsolver.c kernel_fun.c

 Synopsis:
  [Alpha,bias,exitflag,kercnt,access,trnerr,t,NA,UB,LB,History] = 
     new_bsvm2_mex(data,labels,ker,arg,C,solver,tmax,tolabs,tolrel,thlb,cache,verb)

 Input:
  data [dim x num_data] Training vectors.
  labels [1 x num_data] Labels.
  ker [string] Kernel identifier.
  arg [1 x nargs] Kernel argument.
  C [1x1] Regularization constant.
  solver [string] Solver; options are 'mdm','imdm','iimdm'
     'keerthi','kowalczyk'.
  tmax [1x1] Maximal number of iterations.
  tolabs [1x1] Absolute tolerance stopping condition.
  tolrel [1x1] Relative tolerance stopping condition.
  thlb [1x1] Threshold on the lower bound.
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
  access [1x1] Number of requested columns of virtual kernel matrix.
  trnerr [1x1] Training error.
  t [1x1] Number of iterations.
  NA [1x1] Number of non-zero alphas returned by the QP solver.
  UB [1x1] Upper bound on the optimal solution.
  LB [1x1] Lower bound on the optimal solution.
  History [2x(t+1)] UB and LB with respect to number of iterations.

  About: Statistical Pattern Recognition Toolbox
  (C) 1999-2004, Written by Vojtech Franc and Vaclav Hlavac
  <a href="http://www.cvut.cz">Czech Technical University Prague</a>
  <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
  <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

 Modifications:
 28-nov-2004, VF
 26-nov-2004, VF
 24-nov-2004, VF
 20-nov-2004, VF
 31-may-2004, VF
 25-jan-2003, VF
 24-jan-2003, VF
 23-jan-2003, VF
-------------------------------------------------------------------- */

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
/*#include <values.h>*/

#include "kernel_fun.h"

#define INTEGER_MAX  1000000000
#define INDEX(ROW,COL,DIM) ((COL*DIM)+ROW)

#define ABS(x) ((x >=0 ) ? (x) : -(x))
#define KDELTA(A,B) (A==B)
#define KDELTA4(A1,A2,A3,A4) ((A1==A2)||(A1==A3)||(A1==A4)||(A2==A3)||(A2==A4)||(A3==A4))


/* Declaration of global variables */

unsigned long access_cnt;
long num_classes;      
long num_virt_data;  /* number of virtual "single-class" examples */
long num_data;       /* number of input training examples */
double kernel_diag;  /* regularization constant */
double *labels;      /* Pointer to labels */

long Cache_Size;     /* number of cached columns (min 1) */

/* cache (FIFO) for columns of the kernel matrix */
long *cache_index;                  /* indices cached of kernel columns */
long first_kernel_inx;              /* index of first inserted column */
double **kernel_columns;            /* pointers at cached columns */

/* cache for three columns of the virtual kernel matrix */
int first_virt_inx;                 /* index of first used column */
double *virt_columns[3];            /* cache for three columns*/

/* ------------------------------------------------------------
  Returns pointer at a-th column of the kernel matrix.
  This function maintains FIFO cache of kernel columns.
------------------------------------------------------------ */
void *get_kernel_col( long a ) 
{
  double *col_ptr;
  long i;
  long inx;

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

  for( i=0; i < num_data; i++ ) {
    col_ptr[i] = kernel(i,a);
  }

  return( col_ptr );
}

/* ------------------------------------------------------------
  Computes index of input example and its class label from 
  index of virtual "single-class" example.
------------------------------------------------------------ */
void get_indices2( long *index, long *class, long i )
{
   *index = i / (num_classes-1);
 
   *class = (i % (num_classes-1))+1;
   if( *class >= labels[ *index ]) (*class)++;

   return;
}

/* ------------------------------------------------------------
  Retures (a,b)-th element of the virtual kernel matrix 
  of size [num_virt_data x num_virt_data]. 
------------------------------------------------------------ */
double kernel_fce( long a, long b )
{
  double value;
  long i1,c1,i2,c2;

  get_indices2( &i1, &c1, a );
  get_indices2( &i2, &c2, b );

  if( KDELTA4(labels[i1],labels[i2],c1,c2) ) {
    value = (+KDELTA(labels[i1],labels[i2]) 
             -KDELTA(labels[i1],c2)
             -KDELTA(labels[i2],c1)
             +KDELTA(c1,c2)
            )*(kernel( i1, i2 )+1);
  }
  else
  {
    value = 0;
  }

  if(a==b) value += kernel_diag; 

  return( value );
}

/* ------------------------------------------------------------
  Returns pointer at the a-th column of the virtual K matrix.

  (note: the b-th column must be preserved in the cache during 
   updating but b is from (a(t-2), a(t-1)) where a=a(t) and
   thus FIFO with three columns does not have to take care od b.)
------------------------------------------------------------ */
void *get_col( long a, long b )
{
  long i;
  long inx;
  long min_usage; 
  double *col_ptr;
  double *ker_ptr;
  double value;
  long i1,c1,i2,c2;

  access_cnt = access_cnt + 1;

  col_ptr = virt_columns[first_virt_inx++];
  if( first_virt_inx >= 3 ) first_virt_inx = 0;

  get_indices2( &i1, &c1, a );
  ker_ptr = (double*) get_kernel_col( i1 );

  for( i=0; i < num_virt_data; i++ ) {
    get_indices2( &i2, &c2, i );

    if( KDELTA4(labels[i1],labels[i2],c1,c2) ) {
      value = (+KDELTA(labels[i1],labels[i2]) 
               -KDELTA(labels[i1],c2)
               -KDELTA(labels[i2],c1)
               +KDELTA(c1,c2)
              )*(ker_ptr[i2]+1);
    }
    else
    {
      value = 0;
    }

    if(a==i) value += kernel_diag; 

    col_ptr[i] = value;
  }
  
  return( col_ptr );
}


/* -------------------------------------------------------------------
 Main MEX function - interface to Matlab.
-------------------------------------------------------------------- */
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[] )
{
  char solver[20];   /* solver identifier */
  int exitflag;      /* output arg */
  int *err_bit;      /* axiliary cache for computation of trn errors*/
  int buf_len;       /* real length of the solver identifier */
  long i,j ;         /* common use loop variables */
  long inx1, inx2;   
  long NA; 
  long tmax;         /* input arg - max number of iteration */ 
  long t;            /* output arg - number of iterations */
  long verb;         /* input argument */
  double thlb;       /* input arg - threshold on lower bound */
  double C;          /* input arg - regularization const */
  double tolrel;     /* input arg */
  double tolabs;     /* input arg */
  double trnerr;     /* output arg */
  double *tmp_ptr;  
  double *tmp_ptr1;
  double *tmp_ptr2; 
  double *vector_c;  /* auxiliary */ 
  double *Alpha;     /* solution vector */ 
  double *History;   /* output arg */
  double *diagK;     /* cache for diagonal of virtual K matrix */

  /*------------------------------------------------------------------- */
  /* Take input arguments                                               */
  /*------------------------------------------------------------------- */

  if( nrhs != 12) mexErrMsgTxt("Incorrect number of input arguments.");

  dataA = mxGetPr(prhs[0]);   /* pointers at data */
  dataB = dataA;
  dim = mxGetM(prhs[0]);      /* data dimension */
  num_data = mxGetN(prhs[0]); /* number of data */
  labels = mxGetPr(prhs[1]);  /* pointer at data labels */

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
  tmax = mxIsInf( mxGetScalar(prhs[6])) ? INTEGER_MAX : (long)mxGetScalar(prhs[6]); 
  tolabs = mxGetScalar(prhs[7]);   /* abs. precision defining stopping cond*/
  tolrel = mxGetScalar(prhs[8]);   /* rel. precision defining stopping cond*/
  /* threshold on lower bound */
  thlb = mxIsInf( mxGetScalar(prhs[9])) ? mxGetInf() : (double)mxGetScalar(prhs[9]); 

  Cache_Size = (long)mxGetScalar(prhs[10]);  /* cache size */
  if( Cache_Size < 1 ) mexErrMsgTxt("Cache must be greater than 1."); 
  if( Cache_Size > num_data ) Cache_Size = num_data; 

  verb = (long)mxGetScalar(prhs[11]);  /* verbosity on/off */

  /*------------------------------------------------------------------- */
  /* Inicialization (caches, etc.)                                      */
  /*------------------------------------------------------------------- */


  /* constant added to diagonal of separable problem */
  if( C!=0 ) kernel_diag = 1/(2*C); else kernel_diag = 0;

  /* num_classes = max( labels ) */
  num_classes = -INTEGER_MAX; 
  for( i = 0; i < num_data; i++ ) { 
     if( labels[i] > num_classes ) num_classes = (long)labels[i]; 
  }

  /* computes number of virtual "single-class" examples */
  num_virt_data = (num_classes-1)*num_data;

  ker_cnt = 0;    /* counter of kernel evaluations */
  access_cnt = 0;  /* counter for access to the kernel matrix */

  /* allocattes and precomputes diagonal of virtual K matrix */
  diagK = mxCalloc(num_virt_data, sizeof(double));
  if( diagK == NULL ) mexErrMsgTxt("Not enough memory.");
  for(i = 0; i < num_virt_data; i++ ) {
    diagK[i] = kernel_fce(i,i);
  }

  /* allocates memory for kernel cache */
  kernel_columns = mxCalloc(Cache_Size, sizeof(double*));
  if( kernel_columns == NULL ) mexErrMsgTxt("Not enough memory.");
  cache_index = mxCalloc(Cache_Size, sizeof(double));
  if( cache_index == NULL ) mexErrMsgTxt("Not enough memory.");

  for(i = 0; i < Cache_Size; i++ ) {
    kernel_columns[i] = mxCalloc(num_data, sizeof(double));
    if(kernel_columns[i] == NULL) mexErrMsgTxt("Not enough memory.");

    cache_index[i] = -2;
  }
  first_kernel_inx = 0;

  /* allocates memory for three virtual kernel matrix columns */
  for(i = 0; i < 3; i++ ) {
    virt_columns[i] = mxCalloc(num_virt_data, sizeof(double));
    if(virt_columns[i] == NULL) mexErrMsgTxt("Not enough memory.");
  }
  first_virt_inx = 0; 

  /* Solution vector */
  Alpha = mxCalloc(num_virt_data, sizeof(double));
  if( Alpha == NULL ) mexErrMsgTxt("Not enough memory.");

  /* Vector c; for this problem set to zero */
  vector_c = mxCalloc(num_virt_data, sizeof(double));
  if( vector_c == NULL ) mexErrMsgTxt("Not enough memory.");
  for(i = 0; i < num_virt_data; i++ ) vector_c[i] = 0;

  /*------------------------------------------------------------------- */
  /* Call QP solver                                                     */
  /*------------------------------------------------------------------- */

  if ( strcmp( solver, "mdm" )==0 ) {  
     exitflag = gmnp_mdm( &get_col, diagK, vector_c, num_virt_data, tmax, 
         tolabs, tolrel, thlb, Alpha, &t, &History, verb );
  } else if ( strcmp( solver, "imdm" )==0 ) {  
     exitflag = gmnp_imdm( &get_col, diagK, vector_c, num_virt_data, tmax, 
         tolabs, tolrel, thlb, Alpha, &t, &History, verb );
  } else if ( strcmp( solver, "iimdm" )==0 ) {  
     exitflag = gmnp_iimdm( &get_col, diagK, vector_c, num_virt_data, tmax, 
         tolabs, tolrel, thlb, Alpha, &t, &History, verb );
  } else if ( strcmp( solver, "keerthi" )==0 ) {  
     exitflag = gmnp_keerthi( &get_col, diagK, vector_c, num_virt_data, tmax, 
         tolabs, tolrel, thlb, Alpha, &t, &History, verb );
  } else if ( strcmp( solver, "kowalczyk" )==0 ) {  
     exitflag = gmnp_kowalczyk( &get_col, diagK, vector_c, num_virt_data, tmax, 
         tolabs, tolrel, thlb, Alpha, &t, &History, verb );
  } else if ( strcmp( solver, "kozinec" )==0 ) {  
     exitflag = gmnp_kozinec( &get_col, diagK, vector_c, num_virt_data, tmax, 
         tolabs, tolrel, thlb, Alpha, &t, &History, verb );
  } else {
     mexErrMsgTxt("Unknown solver identifier.");
  }

  /*------------------------------------------------------------------- */
  /* Generate outputs                                                   */
  /*------------------------------------------------------------------- */

  /* matrix Alpha [num_classes x num_data] */
  plhs[0] = mxCreateDoubleMatrix(num_classes,num_data,mxREAL);
  tmp_ptr1 = mxGetPr(plhs[0]);

  /* bias vector b [num_classes x 1] */
  plhs[1] = mxCreateDoubleMatrix(num_classes,1,mxREAL);
  tmp_ptr2 = mxGetPr(plhs[1]);

  for( i=0; i < num_classes; i++ ) {
    for( j=0; j < num_virt_data; j++ ) {
       get_indices2( &inx1, &inx2, j );

/*       tmp_ptr1[(inx1*num_classes)+i] += */
/*                      Alpha[j]*(KDELTA(labels[inx1],i+1)+KDELTA(i+1,inx2));*/
       tmp_ptr1[(inx1*num_classes)+i] += 
                      Alpha[j]*(KDELTA(labels[inx1],i+1)-KDELTA(i+1,inx2));
       tmp_ptr2[i] += Alpha[j]*(KDELTA(labels[inx1],i+1)-KDELTA(i+1,inx2));
    }
  }

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
  err_bit = mxCalloc(num_data, sizeof(int));
  if( err_bit == NULL ) mexErrMsgTxt("Not enough memory.");
  for( i=0; i < num_classes; i++ ) {
    for( j=0; j < num_virt_data; j++ ) {
       get_indices2( &inx1, &inx2, j );
       if( ABS(Alpha[j]) > 2*C ) err_bit[inx1] = 1; 
    }
  }

  for( trnerr = 0, i = 0; i < num_data; i++ ) trnerr += err_bit[i];

  trnerr = trnerr/num_data;
  plhs[5] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[5])) = trnerr;

  /* t [1x1] */
  plhs[6] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[6])) = (double)t;

  /* NA [1x1] */
  for( NA = 0, j=0; j < num_virt_data; j++ ) {
     if( Alpha[j] != 0 ) NA++; 
  }

  plhs[7] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[7])) = (double)NA;

  /* UB [1x1] */
  plhs[8] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[8])) = History[INDEX(1,t,2)];

  /* LB [1x1] */
  plhs[9] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[9])) = History[INDEX(0,t,2)];

  /* History [2 x (t+1)] */
  plhs[10] = mxCreateDoubleMatrix(2,t+1,mxREAL);
  tmp_ptr = mxGetPr( plhs[10] );
  for( i = 0; i <= t; i++ ) {
     tmp_ptr[INDEX(0,i,2)] = History[INDEX(0,i,2)];
     tmp_ptr[INDEX(1,i,2)] = History[INDEX(1,i,2)];
  }

  /*------------------------------------------------------------------- */
  /* Free used memory                                                   */
  /*------------------------------------------------------------------- */
  mxFree( vector_c );
  mxFree( Alpha );
  mxFree( History );
  mxFree( diagK );
  for(i = 0; i < Cache_Size; i++ ) mxFree(kernel_columns[i]);
  for(i = 0; i < 3; i++ ) mxFree(virt_columns[i]);
  mxFree( kernel_columns );
  mxFree( cache_index );
}

