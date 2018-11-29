/* --------------------------------------------------------------------
 Evaluation of kernel functions: linear, poly, rbf, sigmoidal.

 About: Statistical Pattern Recognition Toolbox
 (C) 1999-2007, Written by Vojtech Franc and Vaclav Hlavac
 <a href="http://www.cvut.cz">Czech Technical University Prague</a>
 <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
 <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

 Modifications:
 05-mar-07, VF, added example on defining a customized kernel. 
 04-may-04, VF
 04-jun-03, VF
 12-nov-01, VF, sigmoid kenel
 22-sep-01, V. Franc, created.
-------------------------------------------------------------------- */

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <string.h>

/* --- Global variables --------------------------------------------- */

double *dataA;      /* pointer at the fist patterns */
double *dataB;      /* pointer at the second patterns */
long dim;           /* dimension of patterns */
int ker;            /* kernel id (0 - linear, 1 - polynomial, 
                       2 - rbf, 3 - sigmoid */
double *arg1;       /* kernel argument */
long ker_cnt;       /* number of cernel evaluations */

/* "custom" can be replaced by user's own kernel-identifier */
char *kernel_name[] = {"linear","poly","rbf","sigmoid","custom"};


/* -------------------------------------------------------------------
 Computes dot product of a-th and b-th vector.
 c = a'*b
------------------------------------------------------------------- */
double dot_prod( long a, long b)
{
   double c = 0;
   long i;
   for( i = 0; i < dim; i++ ) {
      c += *(dataA+(a*dim)+i) * *(dataB+(b*dim)+i);
   }
   return( c );
}

/* -------------------------------------------------------------------
 Computes dot product of subtraction of a-th and b-th vector.
 c = (a-b)'*(a-b)
------------------------------------------------------------------- */
double sub_dot_prod( long a, long b )
{
   double c = 0;
   long i;
   for( i = 0; i < dim; i++ ) {
      c += (*(dataA+(a*dim)+i) - *(dataB+(b*dim)+i))*
           (*(dataA+(a*dim)+i) - *(dataB+(b*dim)+i));
   }
   return( c );
}

/* --------------------------------------------------------------------
 Converts string kernel identifier to int.
-------------------------------------------------------------------- */
int kernel_id( const mxArray *prhs1 )
{
  int num, i, buf_len;
  char *buf;

  if( mxIsChar( prhs1 ) != 1) return( -1 );

  buf_len  = (mxGetM(prhs1) * mxGetN(prhs1)) + 1;
  buf = mxCalloc( buf_len, sizeof( char ));

  mxGetString( prhs1, buf, buf_len );
  
  num = sizeof( kernel_name )/sizeof( char * );

  for( i = 0; i < num; i++ ) {
    if( strcmp( buf, kernel_name[i] )==0 ) return( i );
  }

  return(-1);
}

/* --------------------------------------------------------------------
 Computes kernel function for a-th and b-th.

 The base address for the 1st argument is dataA and for the 2nd
 argument is dataB.
-------------------------------------------------------------------- */
double kernel( long a, long b )
{
   double c = 0;

   ker_cnt++;

   switch( ker ) {
      /* linear kernel */
      case 0:
         c = dot_prod( a, b );
         break;
      /* polynomial kernel */
      case 1:
         c = pow( (dot_prod( a, b) + arg1[1]), arg1[0] );
         break;
      /* radial basis functions kernel*/
      case 2:
         c = exp( -0.5*sub_dot_prod( a, b)/(arg1[0]*arg1[0]) );
         break;
      /* sigmoid */
      case 3:     
         c = tanh( arg1[0]*dot_prod( a,b) + arg1[1] );
         break;
      /* "custom": here comes definition of user's own kernel. 
      Currently, the linear kernel is used. */
      case 4:     
         c = dot_prod( a, b ); 
         break;
      default:
         c = 0;
   }
   return( c );
}

