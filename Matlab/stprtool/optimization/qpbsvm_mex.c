/*---------------------------------------------------------------------------
 qpbsvm_mex.c: Solver for QP task required for learning SVM without bias term.

 Synopsis:
 [x,exitflag,t,access,History] = qpbsvm_mex(H,f,UB,solver,tmax,tolabs,tolrel,tolKKT,verb,x0,Nabla0)

 Compile: 
  mex qpbsvm_mex.c qpbsvmlib.c

 Description:  
  
   min 0.5*x'*H*x + f'*x  subject to  UB >= x(i) >= 0 for all i
  
 Input:
  H [dim x dim] Symmetric positive semi-definite matrix.
  f [dim x 1] Vector.
  UB [1x1] Upper bound on x(i).
  solver [string] 
      'sca' ... Sequential Coordinatewise Algorithm (Generalize Gauss-Seidel); 
      'scas' ... Greedy variant - udpate variable yielding the best improvement.  
  tmax [1x1] Maximal number of iterations.
  tolabs [1x1] Absolute tolerance stopping condition.
  tolrel [1x1] Relative tolerance stopping condition.
  verb [1x1] If 1 then some info about the training is printed.
  
 optional:
  x0 [dim x 1] Initial solution vector.
  mu0 [dim x 1] mu0 = H*x0 + f 

 Output:
  x [dim x 1] Solution vector.
  exitflag [1x1] Indicates which stopping condition took place:
    t >= tmax                   ->  exit_flag = 0  Number of iterations.
    Q_P-Q_D <= tolabs           ->  exit_flag = 1  Abs. tolerance.
    Q_P-Q_D <= Q_P*tolrel       ->  exit_flag = 2  Relative tolerance.
  t [1x1] Number of iterations.
  access [1x1] Access to elements of the matrix H.
  History [2x(t+1)] Q_P and Q_D with respect to number of iterations.

 Modifications:
 20-nov-2006, VF
 01-nov-2006, VF
-------------------------------------------------------------------- */

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#define INDEX(ROW,COL,NUM_ROWS) ((COL*NUM_ROWS)+ROW)
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

 [x,exitflag,t,access,History] = qpbsvm_mex(H,f,UB,tmax,tolabs,tolrel,verb,x0,Nabla0)

-------------------------------------------------------------------- */
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[] )
{
  char solver[20];   /* solver identifier */
  int exitflag;      /* output arg */
  int verb;          /* input argument -- verbosity */
  int buf_len;       /* real length of the solver identifier */
  long i, j;         /* loop variable */
  long tmax;         /* input arg - max number of iteration */ 
  long t;            /* output arg - number of iterations */
  double tolrel;     /* input arg */
  double tolabs;     /* input arg */
  double tolKKT; 
  double *tmp_ptr;  
  double *x;         /* output arg -- solution*/ 
  double *x0;
  double *Nabla;
  double *History;   /* output arg */
  double *diag_H;    /* diagonal of matrix H */
  double *f;         /* vector f */
  double UB;         
  
  /*------------------------------------------------------------------- */
  /* Take input arguments                                               */
  /*------------------------------------------------------------------- */

  if( nrhs != 11) mexErrMsgTxt("Incorrect number of input arguments.");

  /* matrix H */
  matrix_H = mxGetPr(prhs[0]);
  dim = mxGetM(prhs[0]);

  if(dim != mxGetN(prhs[0])) mexErrMsgTxt("Matrix H mast be squared.");
   
  /* vector f */
  f = mxGetPr(prhs[1]);   
  if((MAX(mxGetM(prhs[1]),mxGetN(prhs[1])) != dim) ||
     (MIN(mxGetM(prhs[1]),mxGetN(prhs[1])) != 1))
      mexErrMsgTxt("Vector f is of wrong size.");
  UB = mxGetScalar(prhs[2]);

  /* string identifier of QP solver to be used */
  if( mxIsChar( prhs[3] ) != 1) mexErrMsgTxt("Solver must be a string.");
  buf_len = (mxGetM(prhs[3]) * mxGetN(prhs[3])) + 1;
  buf_len = (buf_len > 20) ? 20 : buf_len;
  mxGetString( prhs[3], solver, buf_len );

  tmax = mxIsInf( mxGetScalar(prhs[4])) ? INT_MAX : (long)mxGetScalar(prhs[4]);

  tolabs = mxGetScalar(prhs[5]);   /* abs. precision defining stopping cond*/
  tolrel = mxGetScalar(prhs[6]);   /* rel. precision defining stopping cond*/
  tolKKT = mxGetScalar(prhs[7]);   /* rel. precision defining stopping cond*/
  verb = (int)(mxGetScalar(prhs[8]));  /* verbosity on/off */

 /* output "solution" vector alpha [dim x 1] */
  plhs[0] = mxCreateDoubleMatrix(dim,1,mxREAL);
  x = mxGetPr(plhs[0]);
  
  x0 = mxGetPr(prhs[9]);
  for(i=0; i < dim; i++) x[i] = x0[i];

  /* Nabla = H*x + f */
  plhs[5] = mxCreateDoubleMatrix(dim,1,mxREAL);
  Nabla = mxGetPr(plhs[5]);
  tmp_ptr = mxGetPr(prhs[10]); 
  for(i=0; i < dim; i++) Nabla[i] = tmp_ptr[i]; 
    
  if( verb > 0 ) {
    mexPrintf("Settings of QP solver\n");
    mexPrintf("nrhs   : %d\n", nrhs);
    mexPrintf("solver : %s\n", solver);
    mexPrintf("UB     : %f\n", UB );
    mexPrintf("tmax   : %d\n", tmax );
    mexPrintf("tolabs : %f\n", tolabs );
    mexPrintf("tolrel : %f\n", tolrel );
    mexPrintf("dim    : %d\n", dim );
    mexPrintf("verb   : %d\n", verb );
  }

  /*------------------------------------------------------------------- */ 
  /* Inicialization                                                     */
  /*------------------------------------------------------------------- */

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

  if ( strcmp( solver, "sca" ) == 0 ) {  
     exitflag = qpbsvm_sca( &get_col, diag_H, f, UB, dim, tmax, 
                tolabs, tolrel, tolKKT, x, Nabla, &t, &History, verb );
  } else if (strcmp( solver, "scas" ) == 0 ) {  
     exitflag = qpbsvm_scas( &get_col, diag_H, f, UB, dim, tmax, 
                tolabs, tolrel, tolKKT, x, Nabla, &t, &History, verb );
  } else if (strcmp( solver, "scamv" ) == 0 ) {  
     exitflag = qpbsvm_scamv( &get_col, diag_H, f, UB, dim, tmax, 
                tolabs, tolrel, tolKKT, x, Nabla, &t, &History, verb );
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

