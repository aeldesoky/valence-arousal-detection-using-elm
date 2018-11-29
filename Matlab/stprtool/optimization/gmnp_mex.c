/*-----------------------------------------------------------------------
 gmnp_mex.c: MEX-file solver for Generalized Minimal Norml Problem

 Synopsis:
 [alpha,exitflag,t,access,History] = gmnp_mex(H,c,solver,tmax,tolabs,tolrel,thlb,verb)

 Compile: 
  mex gmnp_mex.c gmnpsolver.c

 Description:
   The Generalized Minimal Norm Problem to solve reads
  
   min 0.5*alpha'*H*alpha + c'*alpha

   subject to  sum(alpha) = 1,  alpha(i) >= 0
  
 Input:
  H [dim x dim] Symmetric positive definite matrix.
  c [dim x 1] Vector.
  solver [string] GMNP solver: options are 'mdm', 'imdm', 
     'iimdm','kowalczyk','keerthi','kozinec'.
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
 21-nov-2004, VF
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
  
  /*------------------------------------------------------------------- */
  /* Take input arguments                                               */
  /*------------------------------------------------------------------- */

  if( nrhs != 8) mexErrMsgTxt("Incorrect number of input arguments.");

  /* matrix H */
  matrix_H = mxGetPr(prhs[0]);   
  dim = mxGetM(prhs[0]);     

  if(dim != mxGetN(prhs[0])) mexErrMsgTxt("Matrix H mast be squared.");
   
  /* vector c */
  vector_c = mxGetPr(prhs[1]);   
  if((MAX(mxGetM(prhs[1]),mxGetN(prhs[1])) != dim) ||
     (MIN(mxGetM(prhs[1]),mxGetN(prhs[1])) != 1))
      mexErrMsgTxt("Vector is of wrong size.");

  /* string identifier of QP solver to be used */
  if( mxIsChar( prhs[2] ) != 1) mexErrMsgTxt("Solver must be a string.");
  buf_len = (mxGetM(prhs[2]) * mxGetN(prhs[2])) + 1;
  buf_len = (buf_len > 20) ? 20 : buf_len;
  mxGetString( prhs[2], solver, buf_len );

  /* maximal allowed number of iterations */
  tmax = mxIsInf( mxGetScalar(prhs[3])) ? INT_MAX : (long)mxGetScalar(prhs[3]); 
  tolabs = mxGetScalar(prhs[4]);   /* abs. precision defining stopping cond*/
  tolrel = mxGetScalar(prhs[5]);   /* rel. precision defining stopping cond*/

  /* threshold on lower bound */
  thlb = mxIsInf( mxGetScalar(prhs[6])) ? DBL_MAX : (double)mxGetScalar(prhs[6]); 

  verb = (int)mxGetScalar(prhs[7]);  /* verbosity on/off */

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
     exitflag = gmnp_mdm( &get_col, diag_H, vector_c, dim, tmax, 
         tolabs, tolrel, thlb, alpha, &t, &History, verb );
  } else if ( strcmp( solver, "imdm" ) == 0 ) {  
     exitflag = gmnp_imdm( &get_col, diag_H, vector_c, dim, tmax, 
         tolabs, tolrel, thlb, alpha, &t, &History, verb );
  } else if ( strcmp( solver, "iimdm" ) == 0 ) {  
     exitflag = gmnp_iimdm( &get_col, diag_H, vector_c, dim, tmax, 
         tolabs, tolrel, thlb, alpha, &t, &History, verb );
  } else if ( strcmp( solver, "keerthi" ) == 0 ) {  
     exitflag = gmnp_keerthi( &get_col, diag_H, vector_c, dim, tmax, 
         tolabs, tolrel, thlb, alpha, &t, &History, verb );
  } else if ( strcmp( solver, "kowalczyk" ) == 0 ) {  
     exitflag = gmnp_kowalczyk( &get_col, diag_H, vector_c, dim, tmax, 
         tolabs, tolrel, thlb, alpha, &t, &History, verb );
  } else if ( strcmp( solver, "kozinec" ) == 0 ) {  
     exitflag = gmnp_kozinec( &get_col, diag_H, vector_c, dim, tmax, 
         tolabs, tolrel, thlb, alpha, &t, &History, verb );
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

