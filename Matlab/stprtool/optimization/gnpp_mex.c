/*-----------------------------------------------------------------------
 gnpp_mex.c: MEX-file solver for Generalized Nearest Point Problem

 Synopsis:
 [alpha,exitflag,t,access,History] = gnpp_mex(H,c,y,solver,tmax,tolabs,tolrel,thlb,verb)

 Compile: 
  mex gmnp_mex.c gmnpsolver.c

 Description:
   The Generalized Minimal Nearest Point problem to solve reads
  
   min 0.5*alpha'*H*alpha + c'*alpha

  subject to   sum(x(find(y==1))) = 1, 
               sum(x(find(y==2))) = 1, 
               x >= 0.

 Input:
  H [dim x dim] Symmetric positive definite matrix.
  c [dim x 1] Vector.
  y [dim x 1] Vector of labels 1 or 2.
  solver [string] GMNP solver: options are 'mdm', 'imdm'.
  tmax [1x1] Maximal number of iterations.
  tolabs [1x1] Absolute tolerance stopping condition.
  tolrel [1x1] Relative tolerance stopping condition.
  thlb [1x1] Threshold on lower bound.
  verb [1x1] If 1 then some info about the training is printed.

 Output:
  alpha [dim x 1] Solution vector.
  exitflag [1x1] Indicates which stopping condition was used:
    UB-LB <= tolabs           ->  exit_flag = 1   Abs. tolerance.
    UB-LB <= UB*tolrel        ->  exit_flag = 2   Relative tolerance.
    LB > th                   ->  exit_flag = 3   Threshold on LB.
    t >= tmax                 ->  exit_flag = 0   Number of iterations.
  t [1x1] Number of iterations.
  access [1x1] Access to elements of the matrix H.
  History [2x(t+1)] UB and LB with respect to number of iterations.

 Modifications:
 09-sep-2005, VF
-------------------------------------------------------------------- */

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>


#define INDEX(ROW,COL,DIM) ((COL*DIM)+ROW)
#define MIN(A,B) ((A < B) ? A : B)
#define MAX(A,B) ((A > B) ? A : B)

#define MINUS_INF INT_MIN
#define PLUS_INF  INT_MAX

/* ------------------------------------------------------------*/
/* Declaration of global variables                             */
/* ------------------------------------------------------------*/
double *matrix_H;
long dim;
long access;

/* ------------------------------------------------------------
  Returns pointer at the a-th column of the matrix H.
------------------------------------------------------------ */
void *get_col( long a, long b )
{
  access += dim;
  return( &matrix_H[ dim*a ] );
}


/* -------------------------------------------------------------------
 Main MEX function - interface to Matlab.
-------------------------------------------------------------------- */
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[] )
{
  char solver[20];   /* solver identifier */
  int exitflag;      /* output arg */
  int buf_len;       /* real length of the solver identifier */
  int verb;          /* input argument */
  long i,j ;         /* common use loop variables */
  long inx1, inx2;   
  long num_data;     /* number of input training examples */
  long tmax;         /* input arg - max number of iteration */ 
  long t;            /* output arg - number of iterations */
  double tolrel;     /* input arg */
  double tolabs;     /* input arg */
  double thlb;       /* input arg */
  double *tmp_ptr;  
  double *alpha;     /* output arg */ 
  double *History;   /* output arg */
  double *diag_H;    /* diagonal of matrix H */
  double *vector_c;  /* vector c */
  double *vector_y;  
  double aHa11; 
  double aHa22;
  
  /*------------------------------------------------------------------- */
  /* Take input arguments                                               */
  /*------------------------------------------------------------------- */

  if( nrhs != 9) mexErrMsgTxt("Incorrect number of input arguments.");

  /* matrix H */
  matrix_H = mxGetPr(prhs[0]);   
  dim = mxGetM(prhs[0]);     

  if(dim != mxGetN(prhs[0])) mexErrMsgTxt("Matrix H mast be squared.");
   
  /* vector c */
  vector_c = mxGetPr(prhs[1]);   
  if((MAX(mxGetM(prhs[1]),mxGetN(prhs[1])) != dim) ||
     (MIN(mxGetM(prhs[1]),mxGetN(prhs[1])) != 1))
      mexErrMsgTxt("Vector is of wrong size.");

  /* vector y */
  vector_y = mxGetPr(prhs[2]);   
  if((MAX(mxGetM(prhs[2]),mxGetN(prhs[2])) != dim) ||
     (MIN(mxGetM(prhs[2]),mxGetN(prhs[2])) != 1))
      mexErrMsgTxt("Vector is of wrong size.");

  /* string identifier of QP solver to be used */
  if( mxIsChar( prhs[3] ) != 1) mexErrMsgTxt("Solver must be a string.");
  buf_len = (mxGetM(prhs[3]) * mxGetN(prhs[3])) + 1;
  buf_len = (buf_len > 20) ? 20 : buf_len;
  mxGetString( prhs[3], solver, buf_len );

  /* maximal allowed number of iterations */
  tmax = mxIsInf( mxGetScalar(prhs[4])) ? INT_MAX : (long)mxGetScalar(prhs[4]); 
  tolabs = mxGetScalar(prhs[5]);   /* abs. precision defining stopping cond*/
  tolrel = mxGetScalar(prhs[6]);   /* rel. precision defining stopping cond*/

  /* threshold on lower bound */
  thlb = mxIsInf( mxGetScalar(prhs[7])) ? DBL_MAX : (double)mxGetScalar(prhs[7]); 

  verb = (int)mxGetScalar(prhs[8]);  /* verbosity on/off */

  if( verb == 1 ) {
    mexPrintf("Settings of QP solver\n");
    mexPrintf("solver : %s\n", solver );
    mexPrintf("tmax   : %d\n", tmax );
    mexPrintf("tolabs : %f\n", tolabs );
    mexPrintf("tolrel : %f\n", tolrel );
    mexPrintf("dim    : %d\n", dim );
  }

  /*------------------------------------------------------------------- */ 
  /* Inicialization                                                     */
  /*------------------------------------------------------------------- */

  /* output "solution" vector alpha [dim x 1] */
  plhs[0] = mxCreateDoubleMatrix(dim,1,mxREAL);
  alpha = mxGetPr(plhs[0]);

  /* allocattes and precomputes diagonal of virtual K matrix */
  diag_H = mxCalloc(dim, sizeof(double));
  if( diag_H == NULL ) mexErrMsgTxt("Not enough memory.");
  for(i = 0; i < dim; i++ ) {
    diag_H[i] = matrix_H[dim*i+i];
  }

  /* counter of access to matrix H */
  access = dim;
  

  /*------------------------------------------------------------------- */
  /* Call QP solver                                                     */
  /*------------------------------------------------------------------- */

  if ( strcmp( solver, "mdm" ) == 0 ) {  
     exitflag = gnpp_mdm( &get_col, diag_H, vector_c, vector_y, dim, tmax, 
         tolabs, tolrel, thlb, alpha, &t, &aHa11, &aHa22, &History, verb );
  } else if ( strcmp( solver, "imdm" ) == 0 ) {  
     exitflag = gnpp_imdm( &get_col, diag_H, vector_c, vector_y, dim, tmax, 
         tolabs, tolrel, thlb, alpha, &t, &aHa11, &aHa22, &History, verb );
  } else {
     mexErrMsgTxt("Unknown QP solver identifier!");
  }

  /*------------------------------------------------------------------- */
  /* Generate outputs                                                   */
  /*------------------------------------------------------------------- */

  /* exitflag [1x1] */
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[1])) = (double)exitflag;

  /* t [1x1] */
  plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[2])) = (double)t;

  /* access [1x1] */
  plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[3])) = (double)access;

  /* History [2 x (t+1)] */
  plhs[4] = mxCreateDoubleMatrix(2,t+1,mxREAL);
  tmp_ptr = mxGetPr( plhs[4] );
  for( i = 0; i <= t; i++ ) {
     tmp_ptr[INDEX(0,i,2)] = History[INDEX(0,i,2)];
     tmp_ptr[INDEX(1,i,2)] = History[INDEX(1,i,2)];
  }

  /*------------------------------------------------------------------- */
  /* Free used memory                                                   */
  /*------------------------------------------------------------------- */
  mxFree( History );
  mxFree( diag_H );

}

