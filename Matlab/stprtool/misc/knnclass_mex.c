/*---------------------------------------------------------------------------
 knnclass_mex.c: MEX-file code K-NN classifier.                    

 Compile:  knnclass_mex.c
                                                       
 Synopsis:  

   [tst_labels,dfce] = knnclass_mex(tst_data,trn_data,trn_labels, k)

 Input:
  tst_data [dim x n_tst] data to be classified.
  trn_data [dim x n_trn] training data.
  trn_labels [ 1 x n_trn] labels of training data.
  k [int] number of neighbours.

 Output:
  tst_labels [1 x n_tst] estimated labels of testing data.

 About: (c) Statistical Pattern Recognition Toolbox, (C) 1999-2003,
 Written by Vojtech Franc and Vaclav Hlavac,
 <a href="http://www.cvut.cz">Czech Technical University Prague</a>,
 <a href="http://www.feld.cvut.cz">Faculty of Electrical engineering</a>,
 <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

 Modifications:
 9-may-2003, VF
 18-sep-2002, V.Franc
-------------------------------------------------------------------- */
 
#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#define MAX_INF INT_MAX

#define MAX(A,B)   (((A) > (B)) ? (A) : (B) )
#define MIN(A,B)   (((A) < (B)) ? (A) : (B) ) 
 

/* ==============================================================
 Main MEX function - interface to Matlab.
============================================================== */
void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[] )
{
  long n_tst, n_trn, k, i, j, l, inx, dim;
  double *tst_data, *trn_data;
  double *tst_labels, *trn_labels;
  double *dist, adist, max_dist;
  double a, b;
  long max_inx, best_count, best_label, count, *max_labels;

  /* -- gets input arguments --------------------------- */

  tst_data = mxGetPr( prhs[0]);
  trn_data = mxGetPr(prhs[1]);
  dim = mxGetM(prhs[0]);            /* data dimension */
  n_tst = mxGetN(prhs[0]);              /* number of data */
  n_trn = mxGetN(prhs[1]);              /* number of data */
  trn_labels = mxGetPr(prhs[2]);

  if( dim != mxGetM( prhs[1] )) {
     mexErrMsgTxt("Dimension of training and testing data differs.");
  }

  if( nrhs != 4 ) {
     mexErrMsgTxt("Incorrect number of input arguments.");
  }
  
  k = mxGetScalar(prhs[3]);    

  /*  output labels*/
  plhs[0] = mxCreateDoubleMatrix(1,n_tst,mxREAL);
  tst_labels = mxGetPr(plhs[0] );

  /*--------------------------*/

  if( (dist = mxCalloc(n_trn, sizeof(double))) == NULL) {
      mexErrMsgTxt("Not enough memory for error cache.");
   }

  if( (max_labels = (long*) mxCalloc(k, sizeof(long))) == NULL) {
      mexErrMsgTxt("Not enough memory for error cache.");
   }


  for( i=0; i < n_tst; i++ ) {

    for( j=0; j < n_trn; j++ ) {

      adist = 0;
      for( l=0; l < dim; l++ ) {
        a = *(tst_data+(i*dim)+l);
        b = *(trn_data+(j*dim)+l);
        adist += a*a - 2*a*b + b*b; 
      }

      dist[j] = sqrt(adist);
    }

    for( l=0; l < k; l++) {

      max_dist=MAX_INF;
      for( j=0; j < n_trn; j++ ) {
        if( max_dist > dist[j] ) {
          max_inx = j;
          max_dist = dist[j];
        }
      }
      dist[ max_inx ] = MAX_INF;
      max_labels[l] = trn_labels[max_inx];
    }

    best_count=0;
    for( l=0; l < k; l++) {
      count = 0;
      for( j=0; j < k; j++) {
        if( max_labels[l] == max_labels[j] ) count++;
      }
      if( count > best_count ) {
        best_count = count;
        best_label = max_labels[l]; 
      }
    }    

    tst_labels[i] = best_label;

  }
  
  mxFree( dist );  /* free mem*/
  mxFree( max_labels );  /* free mem*/
}

