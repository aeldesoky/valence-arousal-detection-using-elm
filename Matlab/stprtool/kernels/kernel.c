/* --------------------------------------------------------------------

 kernel.c: MEX-file code for evaluation of kernel functions.

 Compile:  mex kernel.c kernel_fun.c

 Synopsis:
 
  K = kernel( data, ker, arg )

    data [dim x n1] ... Input vectors.
    ker [string] ... Kernel identifier (see kernel_fun.c)
    arg [1 x nargarg] ... Kernel argument(s).

    K [n1 x n1] ... Kernel matrix K[i,j] = kernel(dataA(:,i),dataA(:,j));


  K = kernel( dataA, dataB, ker, arg )

    dataA [dim x n1] ... Matrix A.
    dataB [dim x n2] ... Matrix B.
    ker [string] ... Kernel identifier (see kernel_fun.c)
    arg [1 x nargarg] ... Kernel argument(s).

    K [n1 x n2] ... Kernel matrix K[i,j] = kernel(dataA(:,i),dataB(:,j));


 About: Statistical Pattern Recognition Toolbox
 (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
 <a href="http://www.cvut.cz">Czech Technical University Prague</a>
 <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
 <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

 Modifications:
 4-may-2004, VF
 21-jan-2002, VF
 13-sep-2002, VF
 21-October-2001, V.Franc.
 30-September-2001, V.Franc, created.
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
   long i, j, n1, n2;
   double tmp;
   double *K;
  

   /* K = kernel( data, ker, arg ) */
   /* ------------------------------------------- */
   if( nrhs == 3) 
   {
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

     /* get pointer at input vectors */
     dataA = mxGetPr(prhs[0]);   
     dataB = dataA;
     dim = mxGetM(prhs[0]);      
     n1 = mxGetN(prhs[0]);       

     /* creates output kernel matrix. */
     plhs[0] = mxCreateDoubleMatrix(n1,n1,mxREAL);
     K = mxGetPr(plhs[0]);

     /* computes kernel matrix. */
     for( i = 0; i < n1; i++ ) {
        for( j = i; j < n1; j++ ) {
           tmp = kernel( i, j );
           K[i*n1+j] = tmp; 
           K[j*n1+i] = tmp; /* kernel is symetric */
        }
     }
   } 
   /* K = kernel( dataA, dataB, ker, arg ) */
   /* ------------------------------------------- */
   else if( nrhs == 4)
   {
      /* data matrix [dim x n1 ] */
      if( !mxIsNumeric(prhs[0]) || !mxIsDouble(prhs[0]) ||
        mxIsEmpty(prhs[0])    || mxIsComplex(prhs[0]) )
        mexErrMsgTxt("Input dataA must be a real matrix.");

      /* data matrix [dim x n2 ] */
      if( !mxIsNumeric(prhs[1]) || !mxIsDouble(prhs[1]) ||
        mxIsEmpty(prhs[1])    || mxIsComplex(prhs[1]) )
        mexErrMsgTxt("Input dataB must be a real matrix.");

      /* kernel identifier */
      ker = kernel_id( prhs[2] );
      if( ker == -1 ) 
        mexErrMsgTxt("Improper kernel identifier.");

     /*  get pointer to arguments  */
     arg1 = mxGetPr(prhs[3]);

     /* pointer at patterns */
     dataA = mxGetPr(prhs[0]);    
     dataB = mxGetPr(prhs[1]);    
     dim = mxGetM(prhs[0]);       
     n1 = mxGetN(prhs[0]);        
     n2 = mxGetN(prhs[1]);        

     /* creates output kernel matrix. */
     plhs[0] = mxCreateDoubleMatrix(n1,n2,mxREAL);
     K = mxGetPr(plhs[0]);

     /* computes kenrel matrix. */
     for( i = 0; i < n1; i++ ) {
        for( j = 0; j < n2; j++ ) {
           K[j*n1+i] = kernel( i, j );
        }
     }
   }
   else
   {
      mexErrMsgTxt("Wrong number of input arguments.");
   }

   return;
}
