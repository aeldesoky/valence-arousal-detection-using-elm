/*--------------------------------------------------------------------------
 gsmo_mex.c: Matlab MEX interface for the Generalized SMO solver.

 Synopsis:
   [x,exitflag,t,access,Nabla] = 
          gsmo_mex(H,f,a,b,LB,UB,x0,Nabla0,tmax,tolKKT,verb)

 Compile: 
  mex gsmo_mex.c gsmolib.c

 Description:  
   min Q_P(x) = 0.5*x'*H*x + f'*x  
    x                                      

   s.t.    a'*x = b 
           LB(i) <= x(i) <= UB(i)   for all i=1:n
    
 Input:
  H [n x n] Symmetric positive semidefinite matrix.
  f [n x 1] Vector.
  a [n x 1] Vector which must not contain zero entries.
  b [1 x 1] Scalar.
  LB [n x 1] Lower bound; -inf is allowed.
  UB [n x 1] Upper bound; inf is allowed.
  x0 [n x 1] Initial solution (default zeros).
  Nabla0 [n x 1] Nabla0 = H*x0 + f.
  tolKKT [1 x 1] Determines relaxed KKT conditions (tau in Keerthi's paper).
  verb [1 x 1] if > 0 then prints info every verb-th iterations.
  tmax [1 x 1] Maximal number of iterations.
  
 Output:
  x [n x 1] Solution vector.
  exitflag [1 x 1] Indicates which stopping condition was used:
      relaxed KKT conditions satisfied  ->  exitflag = 1  
      t >= tmax                         ->  exitflag = 0
  t [1x1] Number of iterations.
  access [1x1] Access to entries of the matrix H.
  Nabla [dim x 1] Nabla = H*x+f.

 Modifications:
 29-nov-2006, VF
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

  [x,exitflag,t,access,Nabla] = 
         gsmo_mex(H,f,a,b,LB,UB,x0,Nabla0,tmax,tolKKT,verb);

-------------------------------------------------------------------- */
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[] )
{
  int exitflag;      /* output arg */
  int verb;          /* input argument -- verbosity */
  long i;            /* loop variable */
  long tmax;         /* input arg - max number of iteration */ 
  long t;            /* output arg - number of iterations */
  double tolKKT; 
  double *x;         /* output arg -- solution*/ 
  double *x0;         
  double *Nabla;
  double *Nabla0;
  double *diag_H;    /* diagonal of matrix H */
  double *f;         /* vector f */
  double *a;
  double b;
  double *LB;
  double *UB;
  double fval;
  
  /*------------------------------------------------------------------- */
  /* Take input arguments                                               */
  /*------------------------------------------------------------------- */

  if( nrhs != 11) mexErrMsgTxt("Incorrect number of input arguments.");

  matrix_H = mxGetPr(prhs[0]);
  dim = mxGetM(prhs[0]);
  f = mxGetPr(prhs[1]);   
  a = mxGetPr(prhs[2]);
  b = mxGetScalar(prhs[3]);
  LB = mxGetPr(prhs[4]);
  UB = mxGetPr(prhs[5]);
  x0 = mxGetPr(prhs[6]);
  Nabla0 = mxGetPr(prhs[7]); 
  tmax = mxIsInf( mxGetScalar(prhs[8])) ? INT_MAX : (long)mxGetScalar(prhs[8]);
  tolKKT = mxGetScalar(prhs[9]); 
  verb = (int)(mxGetScalar(prhs[10])); 
    
  if( verb > 0 ) {
    mexPrintf("Settings of QP solver\n");
    mexPrintf("nrhs   : %d\n", nrhs);
    mexPrintf("nlhs   : %d\n", nlhs);
    mexPrintf("tmax   : %d\n", tmax );
    mexPrintf("tolKKT : %f\n", tolKKT );
    mexPrintf("n      : %d\n", dim );
    mexPrintf("verb   : %d\n", verb );
  }

  /*------------------------------------------------------------------- */ 
  /* Inicialization                                                     */
  /*------------------------------------------------------------------- */

  plhs[0] = mxCreateDoubleMatrix(dim,1,mxREAL);
  x = mxGetPr(plhs[0]);
  for(i=0; i < dim; i++) x[i] = x0[i];

  plhs[4] = mxCreateDoubleMatrix(dim,1,mxREAL);
  Nabla = mxGetPr(plhs[4]);
  for(i=0; i < dim; i++) Nabla[i] = Nabla0[i]; 

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

  exitflag = gsmo_algo( &get_col, diag_H, f, a, b, LB, UB, x, Nabla, 
                        dim, tmax, tolKKT, verb, &t );

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

  /*------------------------------------------------------------------- */
  /* Free used memory                                                   */
  /*------------------------------------------------------------------- */
  mxFree( diag_H );
}

/* EOF */
