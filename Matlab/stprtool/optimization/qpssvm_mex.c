/* -----------------------------------------------------------------------
qpssvm_mex.c: MEX-file implementation of QP solver for StructSVM learning
                                                                     
Synopsis:                                                                     
 [x,exitflag,t,access,History] =                                     
            qpssvm_mex(H,f,b,I,x0,tmax,tolabs,tolrel,verb)           

Compile:
 mex qpssvm_mex.c qpssvmlib.c

Description:
 This MEX-file solves the following QP task:
  
   min 0.5*x'*H*x + f'*x
    x

 subject to 
 
   sum(x(find(I==k))) <= b   for all k=1:max(I)
   x >= 0

 where I is a vector of indeices such that unique(I) = 1:max(I).

Input:
  H [n x n] Symmetric positive semidefinite matrix.
  f [n x 1] Vector.
  b [1 x 1] Scalar.
  I [unit16 n x 1] Vector of indices such that unique(I) = 1:max(I);
  x0 [n x 1] Initial solution vector.
  tmax [1 x 1] Maximal number of iterations.
  tolabs [1 x 1] Absolute tolerance stopping condition. 
  tolrel [1 x 1] Relative tolerance stopping condition. 
  verb [1 x 1] if > 0 then prints info every verb-th iterations.

 Output:
  x [n x 1] Solution vector.
  exitflag [1 x 1] Indicates which stopping condition was used:
    UB-LB <= tolabs           ->  exitflag = 1   Abs. tolerance.
    UB-LB <= UB*tolrel        ->  exitflag = 2   Relative tolerance.
    t >= tmax                 ->  exitflag = 0   Number of iterations.
  t [1x1] Number of iterations.
  access [1x1] Access to elements of the matrix H.
  History [2x(t+1)] UB and LB with respect to number of iterations.
                                                                    
 Modifications:                                                      
 20-feb-2006, VF
 17-feb-2006, VF
-------------------------------------------------------------------------*/



#include "mex.h"
#include "string.h"

#define INDEX(ROW,COL,NUM_ROWS) ((COL)*(NUM_ROWS)+(ROW))
#define MIN(A,B) ((A < B) ? (A) : (B))
#define MAX(A,B) ((A > B) ? (A) : (B))

/* -- Global variables --------------------------------------*/
double *matrix_H; /* pointer to the Hessian matrix [dim x dim] */
long dim;         /* dimension of  H */
long access;      /* number of elements of H requested during 
                     optimization */ 

/* ------------------------------------------------------------
  Returns pointer at i-th column of the Hessian  matrix H.
------------------------------------------------------------ */
void *get_col( long i )
{
  access += dim;
  return( &matrix_H[ dim*i ] );
}

/* -------------------------------------------------------------
  The gateway routine.  
------------------------------------------------------------- */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  int exitflag;      /* status of the solution (output arg) */
  long t;            /* number of iterations (output arg) */
  double *x;         /* solution vector (output arg)*/
  double *History;   /* UB and LB history (output arg) */
  int verb;          /* verbosity (input arg)*/
  long tmax;         /* max number of iteration (input arg)*/
  double tolrel;     /* stopping condition (input arg) */
  double tolabs;     /* stopping condition (input arg) */
  double *x0;        /* initial solution (input arg)*/
  double *f;         /* vector f (input arg) */
  double b;          /* scalar b (input arg) */
  uint16_T *I;       /* vector of uint16_T (input arg) */
  double *diag_H;    /* diagonal of matrix H */
  long i ;           /* loop variable */
  double *tmp_ptr;

  /*------------------------------------------------------------------- */
  /* Take input arguments                                               */
  /*         [...] = qpssvm_mex(H,f,b,I,x0,tmax,tolabs,tolrel,verb)     */
  /*------------------------------------------------------------------- */

  if( nrhs != 9) mexErrMsgTxt("Incorrect number of input arguments.");

  /* matrix H */
  matrix_H = mxGetPr(prhs[0]);
  dim = mxGetM(prhs[0]);
  if(dim != mxGetN(prhs[0])) mexErrMsgTxt("Matrix H mast be squared.");

  /* vector f */
  f = mxGetPr(prhs[1]);
  if((MAX(mxGetM(prhs[1]),mxGetN(prhs[1])) != dim) ||
     (MIN(mxGetM(prhs[1]),mxGetN(prhs[1])) != 1))
      mexErrMsgTxt("Vector f is of wrong size.");

  /* vector b */
  b = mxGetScalar(prhs[2]);

  /* vector I */
  I = (uint16_T*)mxGetPr(prhs[3]);
  if((MAX(mxGetM(prhs[3]),mxGetN(prhs[3])) != dim) ||
     (MIN(mxGetM(prhs[3]),mxGetN(prhs[3])) != 1))
      mexErrMsgTxt("Vector I is of wrong size.");

  /* vector x0 */
  x0 = mxGetPr(prhs[4]);
  if((MAX(mxGetM(prhs[4]),mxGetN(prhs[4])) != dim) ||
     (MIN(mxGetM(prhs[4]),mxGetN(prhs[4])) != 1))
      mexErrMsgTxt("Vector x0 is of wrong size.");

  /* maximal allowed number of iterations */
  tmax = mxIsInf( mxGetScalar(prhs[5])) ? INT_MAX : (long)mxGetScalar(prhs[5]);

  /* abs. precision defining stopping cond*/
  tolabs = mxGetScalar(prhs[6]);   
  
  /* rel. precision defining stopping cond*/
  tolrel = mxGetScalar(prhs[7]);   

  /* verbosity parameter */
  verb = (int)mxGetScalar(prhs[8]);  /* verbosity on/off */

  /* print input setting if required */  
  if( verb > 0 ) {
    mexPrintf("Settings of QP solver:\n");
    mexPrintf("tmax   : %d\n", tmax );
    mexPrintf("tolabs : %f\n", tolabs );
    mexPrintf("tolrel : %f\n", tolrel );
    mexPrintf("dim    : %d\n", dim );
    mexPrintf("b      : %f\n", b );
    mexPrintf("verb   : %d\n", verb );
  }     
  
  /*------------------------------------------------------------------- 
     Inicialization                                                     
   ------------------------------------------------------------------- */

  /* solution vector x [dim x 1] */
  plhs[0] = mxCreateDoubleMatrix(dim,1,mxREAL);
  x = mxGetPr(plhs[0]);
  for(i=0; i < dim; i++ ) {
     x[i] = x0[i];       
  }

  /* make diagonal of the Hessian matrix */
  diag_H = mxCalloc(dim, sizeof(double));
  if( diag_H == NULL ) mexErrMsgTxt("Not enough memory.");
  /* to replace with memcpy(void *dest, const void *src, size_t n); */
  for(i = 0; i < dim; i++ ) {
    diag_H[i] = matrix_H[dim*i+i];
  }

  /* counter of access to matrix H */
  access = dim;
  
  /*------------------------------------------------------------------- 
   Call the QP solver.
   -------------------------------------------------------------------*/
/* exitflag = qpssvm_solver( &get_col, diag_H, f, b, I, x, dim, tmax, 
      tolabs, tolrel, &t, &History, verb );   
 */
 
  exitflag = qpssvm_imdm( &get_col, diag_H, f, b, I, x, dim, tmax, 
         tolabs, tolrel, &t, &History, verb );   

  
  /*------------------------------------------------------------------- 
    Set up output arguments                                                   
         [x,exitflag,t,access,History] = qpssvm_mex(...)
  ------------------------------------------------------------------- */

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

  /*------------------------------------------------------------------- 
     Free used memory                                                   
  ------------------------------------------------------------------- */
  mxFree( History );
  mxFree( diag_H );  
  
  return;
}
 
