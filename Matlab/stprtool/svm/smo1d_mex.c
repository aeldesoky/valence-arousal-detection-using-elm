/* --------------------------------------------------------------------
 smo1d_mex.c: MEX-file code for SMO for 1D linear case.

 Compile:  mex smo1d_mex.c

 Synopsis:
   [Alpha,b,nsv,kercnt,trnerr,margin] = smo1d_mex(X, y, C, eps, tol);

 Input:
  X [1 x num_data] Input numbers.
  y [1 x num_data] Input labels.
  C [1x1] SVM regularization constant.
  eps [1x1] Tolerance of KKT-conditions.
  tol [1x1] Minimal change of variables.

 Output:
  Alpha [nsv x 1] Multipliers.
  b [1x1] Bias.
  nsv [1x1] Number os support vectors.
  kercnt [1x1] Number of used dot product (scalar in this case) 
   evaluations.
  trnerr [1x1] Training classification error.
  margin [1x1] Margin of the found classifier.

 About: Statistical Pattern Recognition Toolbox
 (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
 <a href="http://www.cvut.cz">Czech Technical University Prague</a>
 <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
 <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

 Modifications:
 14-may-2004, VF
 15-July-2003, VF
 -------------------------------------------------------------------- */

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* if RANDOM is defined then a random element is used within
 optimization procedure as originally suggested by Platt. */
#define RANDOM

#define ZERO  1e-12

#define MAX(A,B)   (((A) > (B)) ? (A) : (B) )
#define MIN(A,B)   (((A) < (B)) ? (A) : (B) )


/* --- Global variables ---------------------------------------------- */

long   num_data;           /* number of training patterns */
long   dim;                /* dimension */
long   ker_cnt;            /* number of dot product evaluations */
double C;                  /* trade-off constants */
double tolerance=0.001;    /* tolerance in KKT fulfilment  */
double eps=0.001;          /* minimal Lagrangeian change */
double *data;              /* pointer at patterns */
double *target;            /* pointer at labels */
double w;                  /* normal vector */

double *alpha;             /* Lagrange multipliers */
double *b;                 /* Bias (threshold) */


/* ==============================================================
 Implementation of Sequential Minimal Optimizer (SMO)
============================================================== */

/* --------------------------------------------------------------
 Optimizes objective function for i1-th and i2-th pattern.
-------------------------------------------------------------- */
long takeStep( long i1, long i2 ) {
   double y1, y2, s;
   long i;
   double alpha1, alpha2; 
   double a1, a2;
   double E1, E2, L, H, k11, k22, k12, eta, Lobj, Hobj;
   double c1, c2;
   double t;
   double b1, b2, bnew;
   double delta_b;
   double t1, t2;
  
   if( i1 == i2 ) return( 0 );

   alpha1 = alpha[i1];
   y1 = target[i1];
   E1 = w*data[i1] - *b - y1;

   alpha2 = alpha[i2];
   y2 = target[i2];
   E2 = w*data[i2] -*b - y2;

   s = y1 * y2;
   if(s < 0)
   {
      L = MAX(0, alpha2 - alpha1);
      H = MIN(C, C + alpha2 - alpha1);
   }
   else
   {
     L = MAX(0, alpha2 + alpha1 - C );
     H = MIN(C, alpha2 + alpha1);
   }

   if( L == H ) return( 0 );

   k11 = data[i1]*data[i1];
   k12 = data[i1]*data[i2];
   k22 = data[i2]*data[i2];
   ker_cnt += 3;
   eta = 2 * k12 - k11 - k22;

   if( eta < 0 ) {
      a2 = alpha2 + y2 * (E2 - E1) / eta;
      if( a2 < L )
         a2 = L;
      else if( a2 > H )
         a2 = H;
   }
   else {
      c1 = eta/2;
      c2 = y2 * (E1-E2)- eta * alpha2;
      Lobj = c1 * L * L + c2 * L;
      Hobj = c1 * H * H + c2 * H;

      if( Lobj > Hobj+eps )
         a2 = L;
      else if( Lobj < Hobj-eps )
         a2 = H;
      else
         a2 = alpha2;
   }

   if( fabs(a2-alpha2) < eps*(a2+alpha2+eps )) return( 0 );

   a1 = alpha1 - s * (a2 - alpha2 );
   if( a1 < ZERO ) {
      a2 += s * a1;
      a1 = 0;
   }
   else if( a1 > C-ZERO ) {
      t = a1-C;
      a2 += s * t;
      a1 = C;
   }

   if( a1 > 0 && a1 < C )
      bnew = *b + E1 + y1 * (a1 - alpha1) * k11 + y2 * (a2 - alpha2) * k12;
   else {
      if( a2 > 0 && a2 < C )
         bnew = *b + E2 + y1 *(a1 - alpha1)*k12 + y2*(a2 - alpha2) * k22;
      else {
         b1 = *b + E1 + y1 * (a1 - alpha1) * k11 + y2 * (a2 - alpha2) * k12;
         b2 = *b + E2 + y1 * (a1 - alpha1) * k12 + y2 * (a2 - alpha2) * k22;
         bnew = (b1 + b2) / 2;
      }
   }

   delta_b = bnew - *b;
   *b = bnew;

   t1 = y1 * (a1-alpha1);
   t2 = y2 * (a2-alpha2);

   w = w + data[i1]*y1*(a1-alpha[i1]) + data[i2]*y2*(a2-alpha[i2]);

   alpha[i1] = a1;  
   alpha[i2] = a2;  

   return( 1 );
}


/* --------------------------------------------------------------
 Finds the second Lagrange multiplayer to be optimize.
-------------------------------------------------------------- */
long examineExample( long i1 )
{
   double y1, alpha1, E1, r1;
   double tmax;
   double E2, temp;
   long k, i2;
   long k0;

   y1 = target[i1];
   alpha1 = alpha[i1];

   E1 = w*data[i1] - *b - y1;

   r1 = y1 * E1;
   if(( r1 < -tolerance && alpha1 < C )
      || (r1 > tolerance && alpha1 > 0)) {
    /* Try i2 by three ways; if successful, then immediately return 1; */

      for( i2 = (-1), tmax = 0, k = 0; k < num_data; k++ ) {
         if( alpha[k] > 0 && alpha[k] < C ) {
            E2 = w*data[k] - *b - target[k];

            temp = fabs(E1 - E2);
            if( temp > tmax ) {
               tmax = temp;
               i2 = k;
            }
         }
      }
      if( i2 >= 0 ) {
         if( takeStep(i1,i2) )
            return( 1 );
      }

#ifdef RANDOM
      for( k0 = rand(), k = k0; k < num_data + k0; k++ ) {
         i2 = k % num_data;
#else
      for( k = 0; k < num_data; k++) {
         i2 = k;
#endif
         if( alpha[i2] > 0 && alpha[i2] < C ) {
            if( takeStep(i1,i2) )
               return( 1 );
         }
      }

#ifdef RANDOM
      for( k0 = rand(), k = k0; k < num_data + k0; k++ ) {
         i2 = k % num_data;
#else
      for( k = 0; k < num_data; k++) {
         i2 = k;
#endif
         if( takeStep(i1,i2) )
            return( 1 );
      }

   } /* if( ... ) */

   return( 0 );
}

/* --------------------------------------------------------------
 Main SMO optimization cycle.
-------------------------------------------------------------- */
void runSMO( void )
{
   long numChanged = 0;
   long examineAll = 1;
   long k;

   while( numChanged > 0 || examineAll ) {
      numChanged = 0;

      if( examineAll ) {
         for( k = 0; k < num_data; k++ ) {
            numChanged += examineExample( k );
         }
      }
      else {
         for( k = 0; k < num_data; k++ ) {
            if( alpha[k] != 0 && alpha[k] != C )
               numChanged += examineExample( k );
         }
      }

      if( examineAll == 1 )
         examineAll = 0;
      else if( numChanged == 0 )
         examineAll = 1;
   }
}

/* ==============================================================
 Main MEX function - interface to Matlab.
============================================================== */
void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[] )
{
   long i,j ;
   double *labels12, *nsv, *trn_err, *margin;
   double nerr;

   /* ---- check number of input arguments  ------------- */

   if(nrhs != 5 )
      mexErrMsgTxt("Incorrect number of input arguments.");
   if(nlhs < 2)
      mexErrMsgTxt("Not enough output arguments.");


   /* ---- get input arguments  ----------------------- */
   labels12 = mxGetPr(prhs[1]);    /* labels (1,2) */
   data = mxGetPr(prhs[0]);  /* pointer at data */
   dim = mxGetM(prhs[0]);     /* data dimension */
   num_data = mxGetN(prhs[0]);       /* number of data */
   C = mxGetScalar( prhs[2] );
   eps = mxGetScalar( prhs[3] );
   tolerance = mxGetScalar( prhs[4] );

   /* ---- init variables ------------------------------- */   
   ker_cnt=0;      /* num of dot product evaluations  */

   /* allocate memory for targets (labels) (1,-1) */
   if( (target = mxCalloc(num_data, sizeof(double) )) == NULL) {
      mexErrMsgTxt("Not enough memory.");
   }

   /* transform labels12 (1,2) from to targets (1,-1) */
   for( i = 0; i < num_data; i++ ) {
      target[i] = - labels12[i]*2 + 3;
   }

   /* create output variable for bias */
   plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
   b = mxGetPr(plhs[1]); *b= 0;

   /* create vector for Lagrangeians */
   plhs[0] = mxCreateDoubleMatrix(num_data,1,mxREAL);
   alpha = mxGetPr(plhs[0]);

   /* inicialize alpha  */
   for( i = 0; i < num_data; i++ ) {
     alpha[i] = 0;
   }
   w=0;

   /* ---- run SMO ------------------------------------------- */
   runSMO();


   /* ---- outputs ---------------------------------- */
   if( nlhs >= 3 ) {

      /* count number of support vectors */
      plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
      nsv = mxGetPr(plhs[2]);
      *nsv = 0;

      for( i = 0; i < num_data; i++ ) {
         if( alpha[i] > 0) (*nsv)++; 
      }
   }

   if( nlhs >= 4 ) {
     /* number of used iterations */
     plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
     (*mxGetPr(plhs[3])) = (double)ker_cnt;
   }

   if( nlhs >= 5) {

     /* evaluates classification error on traning patterns */
     plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL);
     trn_err = mxGetPr(plhs[4]);
     *trn_err = 0;

     for( i = 0; i < num_data; i++ ) {
       if( target[i]*(w*data[i]-*b) < 0) (*trn_err)++;
     }

     *trn_err = (*trn_err)/(double)num_data;
   }

   if( nlhs >= 6) {
     
      /* compute margin */
      plhs[5] = mxCreateDoubleMatrix(1,1,mxREAL);
      margin = mxGetPr(plhs[5]);
      *margin = 1/sqrt(w*w);
   }

   /* decision function of type <w,x>+b is used */
   *b = -*b;

   /* ----- free memory --------------------------------------- */
   mxFree( target );
}
