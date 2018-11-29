/*-----------------------------------------------------------------------
gnnlslib.c: Library for solving Generalized Non-negative Least Squares problem.
 
 Generalized Non-negative Least Squares Problem to solve is
  
  min 0.5*x'*H*x + f'*x

  subject to  x(i) >= 0 for all i
  
 H [dim x dim] is symmetric positive semi-definite matrix.
 f [dim x 1] is an arbitrary vector.

 The precision of the found solution is given by
 the parameters tmax, tolabs, tolrel and tolKKT which
 define the stopping conditions:

    t >= tmax                   ->  exit_flag = 0  Number of iterations.
    UB-LB <= tolabs             ->  exit_flag = 1  Abs. tolerance.
    UB-LB <= UB*tolrel          ->  exit_flag = 2  Relative tolerance.
    Relaxed KKT cond. satisfied ->  exit_flag = 3  It means that
       mu(i) >= -tolKKT  for all i &&  
       mu(i) <= tolKKT   for i such that x(i) > 0
       where mu = H*x + f                        

 UB ... Upper bound on the optimal solution.
 LB ... Lower bound on the optimal solution.
 t  ... Number of iterations.
 History ... Value of LB and UB wrt. number of iterations.


 The following algorithms are implemented:
 ..............................................

 - Sequential Coordinate-wise Algorithm 
  (Franc, Hlavac: Sequential Coordinate-wise Algorithm
   for Non-negative Least Squares Problem. Research report. 
   CTU-CMP-2005-06. CTU FEL Prague. )
   exitflag = gnnls_sca( &get_col, diag_H, f, dim, tmax, 
               tolabs, tolrel, tolKKT, x, &t, &History, verb );

 - Sequential Coordinate-wise Algorithm with search for the 
   coordinate which gives the highest improvement.
   exitflag = gnnls_scas( &get_col, diag_H, f, dim, tmax, 
               tolabs, tolrel, tolKKT, x, &t, &History, verb );

 Modifications:
 09-seg-2005, VF
 28-aug-2005, VF
-------------------------------------------------------------------- */

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#define HISTORY_BUF 1000000

#define MINUS_INF INT_MIN
#define PLUS_INF  INT_MAX

#define ABS(A) ((A >= 0) ? A : -(A))
#define MIN(A,B) ((A < B) ? (A) : (B))
#define MAX(A,B) ((A >= B) ? (A) : (B))
#define INDEX(ROW,COL,DIM) ((COL*DIM)+ROW)

/* --------------------------------------------------------------

Usage: exitflag = gnnls_sca( &get_col, diag_H, f, dim, tmax, 
               tolabs, tolrel, tolKKT, x, &t, &History, verb )
-------------------------------------------------------------- */
int gnnls_sca(const void* (*get_col)(long,long),
            double *diag_H,
            double *f,
            long dim, 
            long tmax,
            double tolabs,
            double tolrel,
            double tolKKT,
            double *x,
            long  *ptr_t,
            double **ptr_History,
            long verb)
{
  double *Nabla;
  double *History;
  double *col_H;
  double *tmp_ptr;
  double sumx;
  double x_old;
  double delta_x;
  double xHx;
  double UB;
  double LB;
  double xf;
  double min_Nabla;
  double max_Nabla;
  long History_size;
  long t;
  long i;
  long j;
  int exitflag;

  /* ------------------------------------------------------------ */
  /* Initialization                                               */
  /* ------------------------------------------------------------ */

  Nabla = mxCalloc(dim, sizeof(double));
  if( Nabla == NULL ) mexErrMsgTxt("Not enough memory.");

  History_size = (tmax < HISTORY_BUF ) ? tmax+1 : HISTORY_BUF;
  History = mxCalloc(History_size*2,sizeof(double));
  if( History == NULL ) mexErrMsgTxt("Not enough memory.");

  sumx = 0;
  for( i = 0; i < dim; i++ ) {
    if( diag_H[i] > 0 && f[i]*diag_H[i] < 0) {
      sumx += -f[i]/diag_H[i];
    }

    x[i] = 0;
    Nabla[i] = f[i];
  }

  t = 0;
  exitflag = -1;
  while( exitflag == -1 ) 
  {
    t++;     
 
    for(i = 0; i < dim; i++ ) {
      if( diag_H[i] > 0 ) {
        /* variable update */
        x_old = x[i];
        x[i] = MAX(0, x[i] - Nabla[i]/diag_H[i]);
        
        /* update Nabla */
        delta_x = x[i] - x_old;
        if( delta_x != 0 ) {
          col_H = (double*)get_col(i,-1);
          for(j = 0; j < dim; j++ ) {
            Nabla[j] += col_H[j]*delta_x;
          }
        }
      }
    }

    /* compute upper and lower bounds */
    xHx = 0;
    xf = 0;
    min_Nabla = PLUS_INF;
    max_Nabla = MINUS_INF;
    for(j = 0; j < dim; j++ ) {
      xHx += x[j]*(Nabla[j] - f[j]);
      xf += x[j]*f[j];
      if( min_Nabla > Nabla[j]) min_Nabla = Nabla[j];
      if( x[j] > 0 && max_Nabla < Nabla[j]) max_Nabla = Nabla[j];
    }

    UB = 0.5*xHx + xf;
    LB = -sumx*ABS(min_Nabla) - 0.5*xHx;

    /* stopping conditions */
    if(t >= tmax) exitflag = 0;
    else if(UB-LB <= tolabs) exitflag = 1;
    else if(UB-LB <= ABS(UB)*tolrel) exitflag = 2;
    else if(min_Nabla >= -tolKKT && max_Nabla <= tolKKT) exitflag = 3;

    if( verb > 0 && (t % verb == 0 || t==1)) {
      mexPrintf("%d: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
        t, UB, LB, UB-LB,(UB-LB)/ABS(UB));
    }

    /* Store selected values */
    if( t < History_size ) {
      History[INDEX(0,t,2)] = LB;
      History[INDEX(1,t,2)] = UB;
    }
    else {
      tmp_ptr = mxCalloc((History_size+HISTORY_BUF)*2,sizeof(double));
      if( tmp_ptr == NULL ) mexErrMsgTxt("Not enough memory.");
      for( i = 0; i < History_size; i++ ) {
        tmp_ptr[INDEX(0,i,2)] = History[INDEX(0,i,2)];
        tmp_ptr[INDEX(1,i,2)] = History[INDEX(1,i,2)];
      }
      tmp_ptr[INDEX(0,t,2)] = LB;
      tmp_ptr[INDEX(1,t,2)] = UB;
      
      History_size += HISTORY_BUF;
      mxFree( History );
      History = tmp_ptr;
    }

  }


  (*ptr_t) = t;
  (*ptr_History) = History;

  mxFree( Nabla );

  return( exitflag ); 
}

/* --------------------------------------------------------------

Usage: exitflag = gnnls_scas( &get_col, diag_H, f, dim, tmax, 
               tolabs, tolrel, tolKKT, x, &t, &History, verb )
-------------------------------------------------------------- */
int gnnls_scas(const void* (*get_col)(long,long),
            double *diag_H,
            double *f,
            long dim, 
            long tmax,
            double tolabs,
            double tolrel,
            double tolKKT,
            double *x,
            long  *ptr_t,
            double **ptr_History,
            long verb)
{
  double *Nabla;
  double *History;
  double *col_H;
  double *tmp_ptr;
  double sumx;
  double x_old;
  double x_new;
  double max_up;
  double max_x;
  double delta_x;
  double xHx;
  double UB;
  double LB;
  double up;
  double xf;
  double min_Nabla;
  double max_Nabla;
  long History_size;
  long t;
  long i;
  long j;
  long max_i;
  int exitflag;

  /* ------------------------------------------------------------ */
  /* Initialization                                               */
  /* ------------------------------------------------------------ */

  Nabla = mxCalloc(dim, sizeof(double));
  if( Nabla == NULL ) mexErrMsgTxt("Not enough memory.");

  History_size = (tmax < HISTORY_BUF ) ? tmax+1 : HISTORY_BUF;
  History = mxCalloc(History_size*2,sizeof(double));
  if( History == NULL ) mexErrMsgTxt("Not enough memory.");

  sumx = 0;
  for( i = 0; i < dim; i++ ) {
    if( diag_H[i] > 0 && f[i]*diag_H[i] < 0) {
      sumx += -f[i]/diag_H[i];
    }

    x[i] = 0;
    Nabla[i] = f[i];
  }

  t = 0;
  exitflag = -1;
  while( exitflag == -1 ) 
  {
    t++;     
 
    max_up = MINUS_INF;
    for(i = 0; i < dim; i++ ) {
      if( diag_H[i] > 0 ) {
        /* variable update */
        x_old = x[i];
        x_new = MAX(0, x[i] - Nabla[i]/diag_H[i]);
        
        up = -0.5*diag_H[i]*(x_new*x_new-x_old*x_old) -
          (Nabla[i] - diag_H[i]*x_old)*(x_new - x_old);

        if( up > max_up ) {
          max_i = i;
          max_up = up;
          max_x = x_new;
        }
      }
    }

    x_old = x[max_i];
    x[max_i] = max_x;

    /* update Nabla */
    delta_x = max_x - x_old;
    if( delta_x != 0 ) {
      col_H = (double*)get_col(max_i,-1);
      for(j = 0; j < dim; j++ ) {
        Nabla[j] += col_H[j]*delta_x;
      }
    }

    /* compute upper and lower bounds */
    xHx = 0;
    xf = 0;
    min_Nabla = PLUS_INF;
    max_Nabla = MINUS_INF;
    for(j = 0; j < dim; j++ ) {
      xHx += x[j]*(Nabla[j] - f[j]);
      xf += x[j]*f[j];
      if( min_Nabla > Nabla[j]) min_Nabla = Nabla[j];
      if( x[j] > 0 && max_Nabla < Nabla[j]) max_Nabla = Nabla[j];
    }

    UB = 0.5*xHx + xf;
    LB = -sumx*ABS(min_Nabla) - 0.5*xHx;

    /* stopping conditions */
    if(t >= tmax) exitflag = 0;
    else if(UB-LB <= tolabs) exitflag = 1;
    else if(UB-LB <= ABS(UB)*tolrel) exitflag = 2;
    else if(min_Nabla >= -tolKKT && max_Nabla <= tolKKT) exitflag = 3;

    if( verb > 0 && ((t % verb) == 0 || t==1)) {
      mexPrintf("%d: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
        t, UB, LB, UB-LB,(UB-LB)/ABS(UB));
    }

    /* Store selected values */
    if( t < History_size ) {
      History[INDEX(0,t,2)] = LB;
      History[INDEX(1,t,2)] = UB;
    }
    else {
      tmp_ptr = mxCalloc((History_size+HISTORY_BUF)*2,sizeof(double));
      if( tmp_ptr == NULL ) mexErrMsgTxt("Not enough memory.");
      for( i = 0; i < History_size; i++ ) {
        tmp_ptr[INDEX(0,i,2)] = History[INDEX(0,i,2)];
        tmp_ptr[INDEX(1,i,2)] = History[INDEX(1,i,2)];
      }
      tmp_ptr[INDEX(0,t,2)] = LB;
      tmp_ptr[INDEX(1,t,2)] = UB;
      
      History_size += HISTORY_BUF;
      mxFree( History );
      History = tmp_ptr;
    }

  }


  (*ptr_t) = t;
  (*ptr_History) = History;

  mxFree( Nabla );

  return( exitflag ); 
}

