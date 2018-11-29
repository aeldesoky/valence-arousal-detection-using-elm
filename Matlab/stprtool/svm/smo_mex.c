/* --------------------------------------------------------------------
 smo_mex.c: MEX-file for Sequential Minimal Optimizer.

 Compile:  mex smo_mex.c ../kernels/kernel_fun.c

 Synopsis:

  [Alpha,bias,nsv,kercnt,trnerr,margin] =
       smo_mex(data,labels,ker,arg,C,eps,tol,init_Alpha,init_bias )

  Input: 
   data [dim x num_data ] Training vectors.
   labels [1 x num_data] Labels (1 or 2).
   ker [string] Kernel identifier. 
   arg [1 x nargs] Kernel argument(s).  
   C [1x1] or [2 x 1] or [num_data x 1] Regularization constant.
   eps [1x1] SMO parameter (default 0.001).
   tol [1x1] Tolerance of KKT-conditions (default 0.001).
   init_Alpha [num_data x 1] Initial values of optimized Lagrangeians.
   init_bias [1x1] Initial bias value.

  Output:
   Alpha [num_data x 1] Optimized Lagrangians.
   bias [1x1] Bias.

   nsv [1x1] Number of Support Vectors (number of Alpha > ZERO_LIM).
   kercnt [1x1] Number of kernel evaluations.
   trnerr [1x1] Training classification error.
   margin [1x1] Margin.

 About: Statistical Pattern Recognition Toolbox
 (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
 <a href="http://www.cvut.cz">Czech Technical University Prague</a>
 <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
 <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

 Modifications:
 23-may-2004, VF
 14-January-2003, VF
 23-october-2001, V.Franc
 16-October-2001, V.Franc
 27-september-2001, V.Franc, roundig of a2 in takeStep removed.
 23-September-2001, V.Franc, different trade-off C1 and C2.
 22-September-2001, V.Franc, kernel.c used.
 19-September-2001, V.Franc, computation of nsv and nerr added.
 17-September-2001, V.Franc, created.
 -------------------------------------------------------------------- */

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "kernel_fun.h"

/* if RANDOM is defined then a random element is used within
 optimization procedure as originally suggested. */
/* #define RANDOM*/

#define ZERO_LIM   1e-9     /* patterns with alpha > ZERO_LIM are SV */

#define MAX(A,B)   (((A) > (B)) ? (A) : (B) )
#define MIN(A,B)   (((A) < (B)) ? (A) : (B) )

#define C(arg)   (const_C[arg])


/* --- Global variables ---------------------------------------------- */

unsigned long N = 0;       /* number of training patterns */
double *const_C;           /* trade-off constants */
double tolerance=0.001;    /* tolerance in KKT fulfilment  */
double eps=0.001;          /* minimal Lagrangeian change */
double *data;              /* pointer at patterns */
double *target;            /* pointer at labels */
double *error_cache;       /* error cache */

double *alpha;             /* Lagrange multipliers */
double *b;                 /* Bias (threshold) */


/* ==============================================================
 Implementation of Sequential Minimal Optimizer (SMO)
============================================================== */

/* --------------------------------------------------------------
 Computes value of the learned function for k-th pattern.
-------------------------------------------------------------- */
double learned_func( long k )
{
   double s = 0.;
   long i;
   for( i = 0; i < N; i++ ) {
      if( alpha[i] > 0 )
         s += alpha[i]*target[i]*kernel(i,k);
   }
  s -= *b;
  return( s );
}

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
   if( alpha1 > 0 && alpha1 < C(i1) )
      E1 = error_cache[i1];
   else
      E1 = learned_func(i1) - y1;

   alpha2 = alpha[i2];
   y2 = target[i2];
   if( alpha2 > 0 && alpha2 < C(i2) )
      E2 = error_cache[i2];
   else
      E2 = learned_func(i2) - y2;

   s = y1 * y2;
   if(s < 0)
   {
      L = MAX(0, alpha2 - alpha1);
      H = MIN(C(i2), C(i1) + alpha2 - alpha1);
   }
   else
   {
     L = MAX(0, alpha2 + alpha1 - C(i1) );
     H = MIN(C(i2), alpha2 + alpha1);
   }

   if( L == H ) return( 0 );

   k11 = kernel(i1,i1);
   k12 = kernel(i1,i2);
   k22 = kernel(i2,i2);
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
   if( a1 < 0 ) {
      a2 += s * a1;
      a1 = 0;
   }
   else if( a1 > C(i1) ) {
      t = a1-C(i1);
      a2 += s * t;
      a1 = C(i1);
   }

   if( a1 > 0 && a1 < C(i1) )
      bnew = *b + E1 + y1 * (a1 - alpha1) * k11 + y2 * (a2 - alpha2) * k12;
   else {
      if( a2 > 0 && a2 < C(i2) )
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

   for( i = 0; i < N; i++ ) {
     if (0 < alpha[i] && alpha[i] < C(i)) {
        error_cache[i] +=  t1 * kernel(i1,i) + t2 * kernel(i2,i) - delta_b;
     }
   }

   error_cache[i1] = 0;
   error_cache[i2] = 0;

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

   if( alpha1 > 0 && alpha1 < C(i1) )
      E1 = error_cache[i1];
   else
      E1 = learned_func(i1) - y1;

   r1 = y1 * E1;
   if(( r1 < -tolerance && alpha1 < C(i1) )
      || (r1 > tolerance && alpha1 > 0)) {
    /* Try i2 by three ways; if successful, then immediately return 1; */

      for( i2 = (-1), tmax = 0, k = 0; k < N; k++ ) {
         if( alpha[k] > 0 && alpha[k] < C(k) ) {
            E2 = error_cache[k];
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
      for( k0 = rand(), k = k0; k < N + k0; k++ ) {
         i2 = k % N;
#else
      for( k = 0; k < N; k++) {
         i2 = k;
#endif
         if( alpha[i2] > 0 && alpha[i2] < C(i2) ) {
            if( takeStep(i1,i2) )
               return( 1 );
         }
      }

#ifdef RANDOM
      for( k0 = rand(), k = k0; k < N + k0; k++ ) {
         i2 = k % N;
#else
      for( k = 0; k < N; k++) {
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
         for( k = 0; k < N; k++ ) {
            numChanged += examineExample( k );
         }
      }
      else {
         for( k = 0; k < N; k++ ) {
            if( alpha[k] != 0 && alpha[k] != C(k) )
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
   double *labels12, *initAlpha, *nsv, *tmp, *trn_err, *margin;
   double nerr;
   double C1, C2;


   /* ---- get input arguments  ----------------------- */
   if(nrhs < 5)
      mexErrMsgTxt("Not enough input arguments.");

   /* data matrix [dim x N ] */
   if( !mxIsNumeric(prhs[0]) || !mxIsDouble(prhs[0]) ||
       mxIsEmpty(prhs[0])    || mxIsComplex(prhs[0]) )
      mexErrMsgTxt("Input X must be a real matrix.");

   /* vector of labels (1,2) */
   if( !mxIsNumeric(prhs[1]) || !mxIsDouble(prhs[1]) ||
       mxIsEmpty(prhs[1])    || mxIsComplex(prhs[1]) ||
       (mxGetN(prhs[1]) != 1 && mxGetM(prhs[1]) != 1))
      mexErrMsgTxt("Input I must be a real vector.");

   labels12 = mxGetPr(prhs[1]);    /* labels (1,2) */
   dataA = mxGetPr(prhs[0]);  /* pointer at patterns */
   dataB = dataA;
   dim = mxGetM(prhs[0]);     /* data dimension */
   N = mxGetN(prhs[0]);       /* number of data */

   /* kernel identifier */
   ker = kernel_id( prhs[2] );
   if( ker == -1 ) 
     mexErrMsgTxt("Improper kernel identifier.");
      
   /*  get pointer to arguments  */
   arg1 = mxGetPr(prhs[3]);

   /*  one or two real trade-off constant(s)  */
   if( !mxIsNumeric(prhs[4]) || !mxIsDouble(prhs[4]) ||
       mxIsEmpty(prhs[4])    || mxIsComplex(prhs[4]) ||
       (mxGetN(prhs[4]) != 1  && mxGetM(prhs[4]) != 1 ))
      mexErrMsgTxt("Improper input argument C.");
   else {
      /* allocate memory for constant C */
      if( (const_C = mxCalloc(N, sizeof(double) )) == NULL) {
        mexErrMsgTxt("Not enough memory.");
      }

      if( MAX( mxGetN(prhs[4]), mxGetM(prhs[4])) == 1 ) {
        C1 = mxGetScalar(prhs[4]);
        for( i=0; i < N; i++ ) const_C[i] = C1; 
      } else
      if( MAX( mxGetN(prhs[4]), mxGetM(prhs[4])) == 2 ) {
         tmp = mxGetPr(prhs[4]);
         C1 = tmp[0];
         C2 = tmp[1];
         for( i=0; i < N; i++ ) {
           if( labels12[i]==1) const_C[i] = C1; else const_C[i] = C2;
         }
      } else
      if( MAX( mxGetN(prhs[4]), mxGetM(prhs[4])) == N ) { 
         tmp = mxGetPr(prhs[4]);
         for( i=0; i < N; i++ ) const_C[i] = tmp[i]; 
      } else {
        mexErrMsgTxt("Improper argument C.");
      }
   }

   /* real parameter eps */
   if( nrhs >= 6 ) {
      if( !mxIsNumeric(prhs[5]) || !mxIsDouble(prhs[5]) ||
         mxIsEmpty(prhs[5])    || mxIsComplex(prhs[5]) ||
         mxGetN(prhs[5]) != 1  || mxGetM(prhs[5]) != 1 )
         mexErrMsgTxt("Input eps must be a scalar.");
      else
         eps = mxGetScalar(prhs[5]);   /* take eps argument */
   }

   /* real parameter tol */
   if(nrhs >= 7) {
      if( !mxIsNumeric(prhs[6]) || !mxIsDouble(prhs[6]) ||
         mxIsEmpty(prhs[6])    || mxIsComplex(prhs[6]) ||
         mxGetN(prhs[6]) != 1  || mxGetM(prhs[6]) != 1 )
         mexErrMsgTxt("Input tol must be a scalar.");
      else
         tolerance = mxGetScalar(prhs[6]);  /* take tolerance argument */
   }

   /* real vector of Lagrangeian multipliers */
   if(nrhs >= 8) {
      if( !mxIsNumeric(prhs[7]) || !mxIsDouble(prhs[7]) ||
          mxIsEmpty(prhs[7])    || mxIsComplex(prhs[7]) ||
          (mxGetN(prhs[7]) != 1  && mxGetM(prhs[7]) != 1 ))
          mexErrMsgTxt("Input Alpha must be a vector.");
   }

   /* real scalar - bias */
   if( nrhs >= 9 ) {
      if( !mxIsNumeric(prhs[8]) || !mxIsDouble(prhs[8]) ||
         mxIsEmpty(prhs[8])    || mxIsComplex(prhs[8]) ||
         mxGetN(prhs[8]) != 1  || mxGetM(prhs[8]) != 1 )
         mexErrMsgTxt("Input bias must be a scalar.");
   }

   /* ---- init variables ------------------------------- */
   
   ker_cnt = 0;

   /* allocate memory for targets (labels) (1,-1) */
   if( (target = mxCalloc(N, sizeof(double) )) == NULL) {
      mexErrMsgTxt("Not enough memory.");
   }

   /* transform labels12 (1,2) from to targets (1,-1) */
   for( i = 0; i < N; i++ ) {
      target[i] = - labels12[i]*2 + 3;
   }

   /* create output variable for bias */
   plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
   b = mxGetPr(plhs[1]);

   /* take init value of bias if given */
   if( nrhs >= 9 ) {
      *b = -mxGetScalar(prhs[8]);  
   }

   /* allocate memory for error_cache */
   if( (error_cache = mxCalloc(N, sizeof(double) )) == NULL) {
      mexErrMsgTxt("Not enough memory for error cache.");
   }

   /* create vector for Lagrangeians */
   plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
   alpha = mxGetPr(plhs[0]);

   /* if Lagrangeians given then use them as initial values */
   if( nrhs >= 8 ) {
      initAlpha = mxGetPr(prhs[7]);
      for( i = 0; i < N; i++ ) {
         alpha[i] = initAlpha[i];
      }

      /* Init error cache for non-bound multipliers. */
      for( i = 0; i < N; i++ ) {
         if( alpha[i] != 0 && alpha[i] != C(i) ) {
            error_cache[i] = learned_func(i) - target[i];
         }
      }
   }

   /* ---- run SMO ------------------------------------------- */
   runSMO();

   /* ---- outputs  --------------------------------- */
   if( nlhs >= 3 ) {

      /* count number of support vectors */
      plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
      nsv = mxGetPr(plhs[2]);
      *nsv = 0;

      for( i = 0; i < N; i++ ) {
         if( alpha[i] > ZERO_LIM ) (*nsv)++; else alpha[i] = 0;
      }
   }

   if( nlhs >= 4 ) {
     plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
     (*mxGetPr(plhs[3])) = (double)ker_cnt;
   }

   if( nlhs >= 5) {

     /* evaluates classification error on traning patterns */
     plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL);
     trn_err = mxGetPr(plhs[4]);
     nerr = 0;

     for( i = 0; i < N; i++ ) {
        if( target[i] == 1 ) {
           if( learned_func(i) < 0 ) nerr++;
        }
        else
           if( learned_func(i) >= 0 ) nerr++;
     }

     *trn_err = nerr/N;
   }

   if( nlhs >= 6) {
     
      /* compute margin */
      plhs[5] = mxCreateDoubleMatrix(1,1,mxREAL);
      margin = mxGetPr(plhs[5]);
      *margin = 0;
      for( i = 0; i < N; i++ ) {
        for( j = 0; j < N; j++ ) {
           if( alpha[i] > 0 && alpha[j] > 0 )
              *margin += alpha[i]*alpha[j]*target[i]*target[j]*kernel(i,j);
        }
      }

      *margin = 1/sqrt(*margin);
   }

   /* decision function of type <w,x>+b is used */
   *b = -*b;

   /* ----- free memory --------------------------------------- */
   mxFree( error_cache );
   mxFree( target );
}
