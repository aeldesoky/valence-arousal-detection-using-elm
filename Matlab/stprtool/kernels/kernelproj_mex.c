/* --------------------------------------------------------------------
 kernelproj_mex.c: MEX-file code implementing kernel projection.

 Compile:  mex kernelproj_mex.c kernel_fun.c

 Synopsis:
 
  Y = kernelproj_mex(X, Alpha, b, sv_X, ker, arg )

 Input:
   X [dim x num_data] Vectors to be projected.
   Alpha [nsv x new_dim] Multipliers.
   b [new_dim x 1] Bias.
   sv_X [dim x nsv] Support vectors.
   ker [string] Kernel identifier.
   arg [1 x narg] Kernel argument.

 Output:
   Y [new_dim x num_data] Projected data.


 About: Statistical Pattern Recognition Toolbox
 (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
 <a href="http://www.cvut.cz">Czech Technical University Prague</a>
 <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
 <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

 Modifications:
 12-nov-2004, VF
 19-sep-2004, VF
 -------------------------------------------------------------------- */

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include "kernel_fun.h"

/* ==============================================================
 Main MEX function - interface to Matlab.
============================================================== */
void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )
{
   long i, j, k, m;
   long nsv, new_dim, num_data;
   double *Alpha;
   double *b;
   double *Y;
   double k_ij;

  
   ker_cnt = 0;

   /* Y = kernelproj_mex(X, Alpha, b, sv_X, ker, arg ) */
   /* ------------------------------------------- */
   if( nrhs == 6) 
   {
      /* data matrix [dim x num_data] */
      if( !mxIsNumeric(prhs[0]) || !mxIsDouble(prhs[0]) ||
        mxIsEmpty(prhs[0])    || mxIsComplex(prhs[0]) )
        mexErrMsgTxt("Input data must be a real matrix.");

      /* multipliers Alpha [nsv  x new_dim] */
      if( !mxIsNumeric(prhs[1]) || !mxIsDouble(prhs[1]) ||
        mxIsEmpty(prhs[1])    || mxIsComplex(prhs[1]) )
        mexErrMsgTxt("Input Alpha must be a real matrix.");

      /* vector b [nsv  x 1] */
      if( !mxIsNumeric(prhs[2]) || !mxIsDouble(prhs[2]) ||
        mxIsEmpty(prhs[2])    || mxIsComplex(prhs[2]) )
        mexErrMsgTxt("Input b must be a real vector.");

      /* kernel identifier */
      ker = kernel_id( prhs[4] );
      if( ker == -1 ) 
        mexErrMsgTxt("Improper kernel identifier.");
      
     /*  get pointer to arguments  */
     arg1 = mxGetPr(prhs[5]);

     /* get pointer at input vectors */
     dataA = mxGetPr(prhs[0]);   
     Alpha = mxGetPr(prhs[1]);
     b = mxGetPr(prhs[2]);
     dataB = mxGetPr(prhs[3]);  

     /* get data dimensions */ 
     dim = mxGetM(prhs[0]);      
     num_data = mxGetN(prhs[0]);       
     nsv = mxGetM(prhs[1]);
     new_dim = mxGetN(prhs[1]);

     if( mxGetM(prhs[2]) != new_dim)
        mexErrMsgTxt("Number of rows of Alpha must equal to size of vector b.");

     /* creates output kernel matrix. */
     plhs[0] = mxCreateDoubleMatrix(new_dim,num_data,mxREAL);
     Y = mxGetPr(plhs[0]);

     /* computes kernel projection */
     for( i = 0; i < num_data; i++ ) {

       for( k = 0; k < new_dim; k++) { 
         Y[k+i*new_dim] = b[k]; 
       }

       for( j = 0; j < nsv; j++ ) {
         k_ij = kernel(i,j);

         for( k = 0; k < new_dim; k++) { 
           if(Alpha[j+k*nsv] != 0 )
              Y[k+i*new_dim] += k_ij*Alpha[j+k*nsv]; 
         }
       }
     }
   } 
   else
   {
      mexErrMsgTxt("Wrong number of input arguments.");
   }

   return;
}
