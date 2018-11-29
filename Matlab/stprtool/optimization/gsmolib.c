/*-----------------------------------------------------------------------
gsmolib.c: C code implementation of the generalized SMO algorithm.
   
Modifications:
 30-nov-2006, VF
 29-nov-2006, VF
-------------------------------------------------------------------- */

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#define ABS(A) (((A) >= 0) ? (A) : -(A))
#define MIN(A,B) (((A) < (B)) ? (A) : (B))
#define MAX(A,B) (((A) >= (B)) ? (A) : (B))
#define INDEX(ROW,COL,DIM) ((COL*DIM)+ROW)

/* --------------------------------------------------------------

Usage: 
 exitflag = gsmo_algo( &get_col, diag_H, f, a, b, LB, UB, x, Nabla, 
                        dim, tmax, tolKKT, verb, &t );

-------------------------------------------------------------- */
int gsmo_algo(const void* (*get_col)(long,long),
            double *diag_H,
            double *f,
            double *a,
            double b,
            double *LB,
            double *UB,
            double *x,
	        double *Nabla,
            long   dim,
            long   tmax,
            double tolKKT,
            int    verb,
            long   *ptr_t)
{
  double *col_u;
  double *col_v;
  double minF_up;
  double maxF_low;
  double tau;
  double F_i;
  double tau_ub, tau_lb;
  double Q_P;
  long t;
  long i;
  long u, v;
  int exitflag;

  /* ------------------------------------------------------------ */
  /* Initialization                                               */
  /* ------------------------------------------------------------ */

  t = 0;
  exitflag = 0;

  while( exitflag == 0 && t < tmax) 
  {
    t++;     

    /* find the most violating pair of variables */
    minF_up = mxGetInf();
    maxF_low = -mxGetInf();
    for(i = 0; i < dim; i++ ) 
    {
      
      F_i = Nabla[i]/a[i];

      if(LB[i] < x[i] && x[i] < UB[i]) 
      { /* i is from I_0 */
        if( minF_up > F_i) { minF_up = F_i; u = i; }
        if( maxF_low < F_i) { maxF_low = F_i; v = i; }
      } 
      else if((a[i] > 0 && x[i] == LB[i]) || (a[i] < 0 && x[i] == UB[i])) 
      { /* i is from I_1 or I_2 */
        if( minF_up > F_i) { minF_up = F_i; u = i; }
      }
      else if((a[i] > 0 && x[i] == UB[i]) || (a[i] < 0 && x[i] == LB[i])) 
      { /* i is from I_3 or I_4 */
        if( maxF_low < F_i) { maxF_low = F_i; v = i; }
      }
    }

    /* check KKT conditions */
    if( maxF_low - minF_up <= tolKKT )
      exitflag = 1;
    else 
    {
      /* SMO update of the most violating pair */
      col_u = (double*)get_col(u,-1);
      col_v = (double*)get_col(v,-1);

      if( a[u] > 0 ) 
         { tau_lb = (LB[u]-x[u])*a[u]; tau_ub = (UB[u]-x[u])*a[u]; }
      else
         { tau_ub = (LB[u]-x[u])*a[u]; tau_lb = (UB[u]-x[u])*a[u]; }

      if( a[v] > 0 )
         { tau_lb = MAX(tau_lb,(x[v]-UB[v])*a[v]); tau_ub = MIN(tau_ub,(x[v]-LB[v])*a[v]); }
      else
         { tau_lb = MAX(tau_lb,(x[v]-LB[v])*a[v]); tau_ub = MIN(tau_ub,(x[v]-UB[v])*a[v]); }

      tau = (Nabla[v]/a[v]-Nabla[u]/a[u])/
            (diag_H[u]/(a[u]*a[u]) + diag_H[v]/(a[v]*a[v]) - 2*col_u[v]/(a[u]*a[v]));

      tau = MIN(MAX(tau,tau_lb),tau_ub);

      x[u] += tau/a[u];
      x[v] -= tau/a[v];

      /* update Nabla */
      for(i = 0; i < dim; i++ ) 
      {
         Nabla[i] += col_u[i]*tau/a[u] - col_v[i]*tau/a[v];
      } 
    }

    if( verb > 0 && (t % verb == 0 || t==1 || t >= tmax || exitflag > 0)) 
    {
      for(Q_P=0, i = 0; i < dim; i++ ) Q_P += 0.5*(x[i]*Nabla[i]+x[i]*f[i]); /* objective function*/

      mexPrintf("t=%d, KKTviol=%f, tau=%f, tau_lb=%f, tau_ub=%f, Q_P=%f\n",
         t, maxF_low - minF_up, tau, tau_lb, tau_ub, Q_P);
    }

  }  

  (*ptr_t) = t;

  return( exitflag ); 
}



