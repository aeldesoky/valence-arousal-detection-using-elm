/* --------------------------------------------------------------------

 diagker.c: MEX-file code for evaluation of diagonal of kernel matrix.

 Compile:  mex diagker.c kernel_fun.c

 Synopsis:
 
  diagK = diagker( data, ker, arg )

    data [dim x n1] ... Input vectors.
    ker [string] ... Kernel identifier (see kernel_fun.c)
    arg [1 x nargarg] ... Kernel argument(s).

    diagK [n1 x 1] ... Kernel matrix 
       diagK[i] = kernel(dataA(:,i),dataA(:,i));


 About: Statistical Pattern Recognition Toolbox
 (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
 <a href="http://www.cvut.cz">Czech Technical University Prague</a>
 <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
 <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

 Modifications:
 5-may-2004, VF
 20-jun-2003, VF
 -------------------------------------------------------------------- */

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "kernel_fun.h"

/* ==============================================================
 Main MEX function - interface to Matlab.
============================================================== */
void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[] )
{
   long i, num_data;
   double *diagK;

   /* K = diagker( data, ker, arg ) */
   /* ------------------------------------------- */
   if( nrhs != 3) 
     mexErrMsgTxt("Improper number of input arguments.");
     
   /* data matrix [dim x n1] */
   if( !mxIsNumeric(prhs[0]) || !mxIsDouble(prhs[0]) ||
     mxIsEmpty(prhs[0])    || mxIsComplex(prhs[0]) )
     mexErrMsgTxt("Input data must be a real matrix.");

   /* kernel identifier */
   ker = kernel_id( prhs[1] );
   if( ker == -1 ) 
     mexErrMsgTxt("Improper kernel identifier.");
      
   /*  get pointer to arguments  */
   arg1 = mxGetPr(prhs[2]);

   dataA = mxGetPr(prhs[0]);    
   dataB = dataA;
   dim = mxGetM(prhs[0]);       
   num_data = mxGetN(prhs[0]);  

   /* creates output kernel matrix. */
   plhs[0] = mxCreateDoubleMatrix(num_data,1,mxREAL);
   diagK = mxGetPr(plhs[0]);

   /* computes kenrel matrix. */
   for( i = 0; i < num_data; i++ ) {
     diagK[i] = kernel(i,i);
   }
}
