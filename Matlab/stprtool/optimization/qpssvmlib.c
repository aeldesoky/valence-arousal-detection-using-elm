/*-----------------------------------------------------------------------
qpssvmlib.c: Library of solvers for QP task required in StructSVM learning.

Synopsis:

  exitflag = qpssvm_solver( &get_col, diag_H, f, b, I, x, n, tmax, 
             tolabs, tolrel, &t, &History, verb );   

Description:
 
 It solves the following QP task:
  
   min 0.5*x'*H*x + f'*x
    x

 subject to 
 
   sum(x(find(I==k))) <= b   for all k=1:max(I)
   x >= 0

 where I is a set of positive indices from (1 to max(I)).

 A precision of the found solution is given by the parameters tmax, 
 tolabs and tolrel which define the stopping conditions:
 
 UB-LB <= tolabs      ->  exitflag = 1   Abs. tolerance.
 UB-LB <= UB*tolrel   ->  exitflag = 2   Relative tolerance.
 t >= tmax            ->  exitflag = 0   Number of iterations.

 UB ... Upper bound on the optimal solution, i.e., Q_P.
 LB ... Lower bound on the optimal solution, i.e., Q_D.
 t  ... Number of iterations.


Inputs/Outputs:

 const void* (*get_col)(long) retunr pointer to i-th column of H
 diag_H [double n x n] diagonal of H.
 f [double n x 1] is an arbitrary vector.
 b [double 1 x 1] scalar
 I [uint16_T n x 1] Indices (1..max(I)); max(I) <= n
 x [double n x 1] solution vector (inital solution).
 n [long 1 x 1] dimension of H.
 tmax [long 1 x 1] Max number of steps.
 tolrel [double 1 x 1] Relative tolerance.
 tolabs [double 1 x 1] Absolute tolerance.
 t [long 1 x 1] Number of iterations.
 History [double 2 x t] Value of LB and UB wrt. number of iterations.
 verb [int 1 x 1] if > 0 then prints info every verb-th iteation.

 For more info refer to TBA

 Modifications:
 20-Feb-2006, VF
 18-feb-2006, VF

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

#define ABS(A) (((A) >= 0) ? (A) : (-A))
#define MIN(A,B) (((A) < (B)) ? (A) : (B))
#define MAX(A,B) (((A) > (B)) ? (A) : (B))
#define INDEX(ROW,COL,DIM) ((COL*DIM)+ROW)


/* --------------------------------------------------------------
 QPSSVM solver 

 Usage: exitflag = qpssvm_solver( &get_col, diag_H, f, b, I, x, n, tmax, 
         tolabs, tolrel, &t, &History, verb );   
-------------------------------------------------------------- */
int qpssvm_solver(const void* (*get_col)(long),
                  double *diag_H,
                  double *f,
                  double b,
                  uint16_T *I,
                  double *x,
                  long n,
                  long tmax,
                  double tolabs,
                  double tolrel,
                  long *ptr_t,
                  double **ptr_History,
                  long verb)
{
  double *Hx;
  double *d;
  double *History;
  double *col_u, *col_v;
  double *tmp_ptr;
  double LB;
  double UB;
  double tmp;
  double yu;
  double den1;
  double num1;
  double tau1;
  double improv1;
  double den2;
  double num2;
  double tau2;
  double improv2;
  double tmp_num;
  double tmp_den;
  double tau;
  double delta;
  double sumx;
  long m;
  long t;
  long u;
  long v;
  long k;
  long i, j;
  long History_size;
  int exitflag;
  
  /* ------------------------------------------------------------ 
    Initialization                                               
  ------------------------------------------------------------ */

  /* count cumber of constraints */
  for( i=0, m=0; i < n; i++ ) m = MAX(m,I[i]);

  /* alloc Hx [n x m] */
  Hx = mxCalloc(m*n, sizeof(double));
  if( Hx == NULL ) mexErrMsgTxt("Not enough memory.");

  /* alloc History [2 x HISTORY_BUF] */
  History_size = (tmax < HISTORY_BUF ) ? tmax+1 : HISTORY_BUF;
  History = mxCalloc(History_size*2,sizeof(double));
  if( History == NULL ) mexErrMsgTxt("Not enough memory.");

  /* alloc d [n x 1] */
  d = mxCalloc(n, sizeof(double));
  if( d == NULL ) mexErrMsgTxt("Not enough memory.");
 
  /* Hx = zeros(n,m);
  for k=1:m,
    inx = find(I==k);
    Hx(:,k) = H(:,inx)*x(inx);
  end
  */
  for( i=0; i < n; i++ ) {
    if( x[i] > 0 ) {
      u = I[i]-1;
      col_u = (double*)get_col(i);
      for( j=0; j < n; j++ ) {
        Hx[INDEX(j,u,n)] += col_u[j]*x[i];      
      }
    }
  }
  
  /* d = sum(Hx,2) + f; */
  for( i=0; i < n; i++ ) {
    for( j=0; j < m; j++ ) {
      d[i] += Hx[INDEX(i,j,n)]; 
    }
    d[i] += f[i];
  }

  /* UB = 0.5*x'*(f+d); */
  /* LB = 0.5*x'*(f-d); */
  for( i=0, UB = 0, LB=0; i < n; i++) {
    UB += x[i]*(f[i]+d[i]);
    LB += x[i]*(f[i]-d[i]);
  }
  UB = 0.5*UB;
  LB = 0.5*LB;

  /*
  for k=1:m,
    tmp = min(d(find(I==k)));
    if tmp < 0, LB = LB + b*tmp; end
  end
  */
  for( i=0; i < m; i++ ) {
    for( j=0, tmp = PLUS_INF; j < n; j++ ) {
      if( I[j]-1 == i ) tmp = MIN(tmp, d[j]);
    }
    if( tmp < 0) LB += b*tmp;
  }

  exitflag = 0;
  t = 0;
  History[INDEX(0,0,2)] = LB;
  History[INDEX(1,0,2)] = UB;


  /* -- Main loop ---------------------------------------- */
  while( (exitflag == 0) && (t < tmax)) 
  {
    t++;

    exitflag = 1;
    for( k=0; k < m; k++ ) 
    {       
      /*
      inx = find(I==k);
      [tmp,u] = min(d(inx)); u = inx(u);
      */
      for( i=0, tmp = PLUS_INF, delta = 0; i < n; i++ ) {
        if( I[i]-1 == k) {
          delta += x[i]*d[i];
          if( tmp > d[i] ) {
            tmp = d[i];
            u = i;
          }
        }
      }

      /* if d(u) < 0, yu = b; else yu = 0; end  */
      if( d[u] < 0) yu = b; else yu = 0;
     
      /* delta = x(inx)'*d(inx) - yu*d(u); */
      delta -= yu*d[u];
      
      if( delta > tolabs/m && delta > tolrel*ABS(UB)/m) 
      {
         exitflag = 0;
         col_u = (double*)get_col(u);
      
         /* -- Kozinec like update ------ */

         /*
         y = x; y(inx) = 0; y(u) = yu;         
         % den1 = (x-y)'*H*(x-y) 
         den1 = (x(inx)-y(inx))'*(Hx(inx,k)-yu*H(inx,u));
         num1 = (x-y)'*d;
         */

         for( i=0, sumx = 0, den1 = 0, num1 = 0; i < n; i++ ) {
           sumx += x[i];
           if( i == u ) {
             num1 += (x[i]-yu)*d[i];
             den1 += (x[i]-yu)*(Hx[INDEX(i,k,n)]-yu*col_u[i]);
           } else if( I[i]-1 == k) {
             num1 += x[i]*d[i];
             den1 += x[i]*(Hx[INDEX(i,k,n)]-yu*col_u[i]);
           } 
         }

         tau1 = MIN(1, num1/den1);
         if( tau1 < 1 ) 
           improv1 = num1*num1/(2*den1); 
         else {
         /* Improv1 = 0.5*x'*(d+f) - 0.5*y'*(yu*H(:,u)+(d-f)-Hx(:,k)) - f'*y; */
           for( i = 0, improv1 = 0; i < n; i++ ) {
             if( i == u ) {
               improv1 += 0.5*x[i]*(d[i]+f[i]) 
                  - 0.5*yu*(yu*col_u[i]+d[i]-f[i]-Hx[INDEX(i,k,n)]) - f[i]*yu;
             } else if( I[i]-1 == k ) {
               improv1 += 0.5*x[i]*(d[i]+f[i]);
             } else {
               improv1 += 0.5*x[i]*(d[i]+f[i]) 
                - 0.5*x[i]*(yu*col_u[i]+d[i]-f[i]-Hx[INDEX(i,k,n)]) - f[i]*x[i];
             }
           }
         }

         /* -- MDM like update --------- */
         improv2 = MINUS_INF;
         if( sumx > 0) {
           for(i = 0; i < n; i++ ) {
             if( (I[i]-1 == k) && (i != u) && (x[i] > 0)) {
                
               tmp_num = d[i] - d[u];
               tmp_den = diag_H[u] - 2*col_u[i] + diag_H[i];
               if( tmp_den > 0 ) {
                 tau = MIN(1,tmp_num/(x[i]*tmp_den));
                 if( tau < 1 ) {
                   tmp = tmp_num*tmp_num/(2*tmp_den);
                 } else {
               /* tmp = x(i)*(d(i)-d(u))-0.5*x(i)^2*(H(u,u)-2*H(i,u)+H(i,i)); */
                   tmp = x[i]*tmp_num-0.5*x[i]*x[i]*tmp_den;
                 }

                 if( tmp > improv2 ) {
                   improv2 = tmp;
                   tau2 = tau;
                   v = i;
                 }
               }
             }
           }
         }

         /* -- Apply the better line segment -------------- */
         if( improv1 > improv2 ) {
           /* 
            d = d + tau1*(yu*H(:,u)-Hx(:,k));
            Hx(:,k) = Hx(:,k)*(1-tau1) + tau1*H(:,u)*yu;

            x(setdif(inx,u)) = x(setdiff(inx,u))*(1-tau1);
            x(u) = x(u)*(1-tau1) + yu*tau1;
           */

           for( i = 0; i < n; i++ ) {
             d[i] += tau1*(yu*col_u[i]-Hx[INDEX(i,k,n)]);
             Hx[INDEX(i,k,n)] = Hx[INDEX(i,k,n)]*(1-tau1) + tau1*yu*col_u[i];

             if( i == u) 
               x[i] = x[i]*(1-tau1) + tau1*yu;
             else if( I[i]-1 == k ) 
               x[i] = x[i]*(1-tau1);
           }
         } 
         else 
         {
           col_v = (double*)get_col(v);
           for( i = 0; i < n; i++ ) {             
             tmp = x[v]*tau2*(col_u[i]-col_v[i]);
             Hx[INDEX(i,k,n)] += tmp;
             d[i] += tmp;
           }           

           x[u] += tau2*x[v];
           x[v] -= tau2*x[v];
         }

         /* mexPrintf("t=%d,k=%d, u=%d, tau1=%f, den1=%f, num1=%f, delta=%f\n", 
             t,k,u,tau1,den1,num1,delta);*/

         /* -- Update the upper bound ---------------------- */
         for( i=0, UB = 0; i < n; i++) {
           UB += x[i]*(f[i]+d[i]);
         }
         UB = 0.5*UB;
      }
    }

    /* -- Computing LB --------------------------------------*/

    /*
    LB = 0.5*x'*(f-d);   % LB = -0.5*x'*H*x;
    for k=1:n,
      tmp = min(d(find(I==k)));
      if tmp < 0, LB = LB + b*tmp; end
    end */
    
    for( i=0, LB=0; i < n; i++) {
       LB += 0.5*x[i]*(f[i]-d[i]);
    }

    for( i=0; i < m; i++ ) { 
      for( j=0, tmp = PLUS_INF; j < n; j++ ) {
        if( I[j]-1 == i ) tmp = MIN(tmp, d[j]);
      }
      if( tmp < 0) LB += b*tmp;
    }

    /* Store LB and UB */
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

    if( verb > 0 && (exitflag > 0 || (t % verb)==0 )) {
       mexPrintf("%d: UB=%.10f, LB=%.10f, UB-LB=%.10f, (UB-LB)/|UB|=%.10f \n",
        t, UB, LB, UB-LB, (UB!=0) ? (UB-LB)/ABS(UB) : 0);      
    }    

  }

  /* -- Find which stopping consition has been used -------- */
  if( UB-LB < tolabs ) exitflag = 1;
  else if(UB-LB < ABS(UB)*tolrel ) exitflag = 2;
  else exitflag = 0;

  /*----------------------------------------------------------   
    Set up outputs                                          
  ---------------------------------------------------------- */
  (*ptr_t) = t;
  (*ptr_History) = History;

  /*---------------------------------------------------------- 
    Clean up
  ---------------------------------------------------------- */
  mxFree( Hx );
  mxFree( d );
  
  return( exitflag ); 

}


/* --------------------------------------------------------------
 QPSSVM solver 

 Usage: exitflag = qpssvm_imdm( &get_col, diag_H, f, b, I, x, n, tmax, 
         tolabs, tolrel, &t, &History, verb );   
-------------------------------------------------------------- */
int qpssvm_imdm(const void* (*get_col)(long),
                  double *diag_H,
                  double *f,
                  double b,
                  uint16_T *I,
                  double *x,
                  long n,
                  long tmax,
                  double tolabs,
                  double tolrel,
                  long *ptr_t,
                  double **ptr_History,
                  long verb)
{
  double *x_nequ;
/*  double *Hx;*/
  double *d;
  double *History;
  double *col_u, *col_v;
  double *tmp_ptr;
  double LB;
  double UB;
  double tmp;
  double improv;
  double tmp_num;
  double tmp_den;
  double tau;
  double delta;
  double yu;
  long *inx;
  long *nk;
  long m;
  long t;
  long u;
  long v;
  long k;
  long i, j;
  long History_size;
  int exitflag;
  
  /* ------------------------------------------------------------ 
    Initialization                                               
  ------------------------------------------------------------ */

  /* count cumber of constraints */
  for( i=0, m=0; i < n; i++ ) m = MAX(m,I[i]);

  /* alloc and initialize x_nequ */
  x_nequ = mxCalloc(m, sizeof(double));
  if( x_nequ == NULL ) mexErrMsgTxt("Not enough memory.");

  /* alloc Inx */
  inx = mxCalloc(m*n, sizeof(double));
  if( inx == NULL ) mexErrMsgTxt("Not enough memory.");

  nk = mxCalloc(m, sizeof(double));
  if( nk == NULL ) mexErrMsgTxt("Not enough memory.");

  for( i=0; i < m; i++ ) x_nequ[i] = b;
  for( i=0; i < n; i++ ) {
     k = I[i]-1;
     x_nequ[k] -= x[i];
     inx[INDEX(nk[k],k,n)] = i;
     nk[k]++;
  }
    
  /* alloc History [2 x HISTORY_BUF] */
  History_size = (tmax < HISTORY_BUF ) ? tmax+1 : HISTORY_BUF;
  History = mxCalloc(History_size*2,sizeof(double));
  if( History == NULL ) mexErrMsgTxt("Not enough memory.");

  /* alloc d [n x 1] */
  d = mxCalloc(n, sizeof(double));
  if( d == NULL ) mexErrMsgTxt("Not enough memory.");
 
  /* d = H*x + f; */
  for( i=0; i < n; i++ ) {
    if( x[i] > 0 ) {
      col_u = (double*)get_col(i);
      for( j=0; j < n; j++ ) {
          d[j] += col_u[j]*x[i];
      }
    }
  }
  for( i=0; i < n; i++ ) d[i] += f[i];
  
  /* UB = 0.5*x'*(f+d); */
  /* LB = 0.5*x'*(f-d); */
  for( i=0, UB = 0, LB=0; i < n; i++) {
    UB += x[i]*(f[i]+d[i]);
    LB += x[i]*(f[i]-d[i]);
  }
  UB = 0.5*UB;
  LB = 0.5*LB;

  /*
  for k=1:m,
    tmp = min(d(find(I==k)));
    if tmp < 0, LB = LB + b*tmp; end
  end
  */
  
  for( i=0; i < m; i++ ) {
    for( j=0, tmp = PLUS_INF; j < nk[i]; j++ ) {
      tmp = MIN(tmp, d[inx[INDEX(j,i,n)]]);
    }
    if( tmp < 0) LB += b*tmp;
  }
  
  /*
  for( i=0; i < m; i++ ) {
    for( j=0, tmp = PLUS_INF; j < n; j++ ) {
      if( I[j]-1 == i ) tmp = MIN(tmp, d[j]);
    }
    if( tmp < 0) LB += b*tmp;
  }*/

  exitflag = 0;
  t = 0;
  History[INDEX(0,0,2)] = LB;
  History[INDEX(1,0,2)] = UB;


  /* -- Main loop ---------------------------------------- */
  while( (exitflag == 0) && (t < tmax)) 
  {
    t++;

    exitflag = 1;
    for( k=0; k < m; k++ ) 
    {       
      /*
      inx = find(I==k);
      [tmp,u] = min(d(inx)); u = inx(u);
      */
        
     for( j=0,tmp = PLUS_INF, delta = 0; j < nk[k]; j++ ) {
        i = inx[INDEX(j,k,n)];
/*      for( i=0, tmp = PLUS_INF, delta = 0; i < n; i++ ) {
        if( I[i]-1 == k) {*/
        delta += x[i]*d[i];
        if( tmp > d[i] ) {
          tmp = d[i];
          u = i;
        }
      }

      /* if d(u) < 0, yu = b; else yu = 0; end  */
      if( d[u] < 0) yu = b; else yu = 0;
     
      /* delta = x(inx)'*d(inx) - yu*d(u); */
      delta -= yu*d[u];
            
      if( delta > tolabs/m && delta > tolrel*ABS(UB)/m) 
      {
         exitflag = 0;
         
         if( yu > 0 ) 
         {
           col_u = (double*)get_col(u);      

           improv = MINUS_INF;
           for( j=0; j < nk[k]; j++ ) {
             i = inx[INDEX(j,k,n)];
           
/*           for(i = 0; i < n; i++ ) {
             if( (I[i]-1 == k) && (i != u) && (x[i] > 0)) {              */
             if(x[i] > 0) {             
               
               tmp_num = x[i]*(d[i] - d[u]); 
               tmp_den = x[i]*x[i]*(diag_H[u] - 2*col_u[i] + diag_H[i]);
               if( tmp_den > 0 ) {
                 if( tmp_num < tmp_den ) {
                    tmp = tmp_num*tmp_num / tmp_den;
                 } else {
                    tmp = tmp_num - 0.5 * tmp_den;
                 }
               }
               if( tmp > improv ) {
                 improv = tmp;
                 tau = MIN(1,tmp_num/tmp_den);
                 v = i;
               }
             }
           }

           tmp_num = -x_nequ[k]*d[u];
           if( tmp_num > 0 ) {
             tmp_den = x_nequ[k]*x_nequ[k]*diag_H[u];
             if( tmp_den > 0 ) {
               if( tmp_num < tmp_den ) {
                 tmp = tmp_num*tmp_num / tmp_den;
               } else {
                   tmp = tmp_num - 0.5 * tmp_den;
               }
             }
           } else {
             tmp = MINUS_INF; 
           }
           
           if( tmp > improv ) {
              tau = MIN(1,tmp_num/tmp_den);
              for( i = 0; i < n; i++ ) {             
                d[i] += x_nequ[k]*tau*col_u[i];
              }
             x[u] += tau*x_nequ[k];
             x_nequ[k] -= tau*x_nequ[k];
               
           } else {
            
             /* updating with the best line segment */
             col_v = (double*)get_col(v);
             for( i = 0; i < n; i++ ) {             
               d[i] += x[v]*tau*(col_u[i]-col_v[i]);
             }

             x[u] += tau*x[v];
             x[v] -= tau*x[v];
           }
         }
         else
         {
           improv = MINUS_INF;
           for( j=0; j < nk[k]; j++ ) {
             i = inx[INDEX(j,k,n)];
           
/*           for(i = 0; i < n; i++ ) {
             if( (I[i]-1 == k) && (x[i] > 0)) {*/
             if( x[i] > 0 && d[i] > 0) {
                
               tmp_num = x[i]*d[i]; 
               tmp_den = x[i]*x[i]*diag_H[i];
               if( tmp_den > 0 ) {
                 if( tmp_num < tmp_den ) {
                    tmp = tmp_num*tmp_num / tmp_den;
                 } else {
                    tmp = tmp_num - 0.5 * tmp_den;
                 }
               }
               if( tmp > improv ) {
                 improv = tmp;
                 tau = MIN(1,tmp_num/tmp_den);
                 v = i;
               }
             }    
           }

           /* updating with the best line segment */
           col_v = (double*)get_col(v);
           for( i = 0; i < n; i++ ) {             
             d[i] -= x[v]*tau*col_v[i];
           }

           x_nequ[k] += tau*x[v];
           x[v] -= tau*x[v];         
         }
                    
/*         for( i=0, UB = 0; i < n; i++) {
            UB += x[i]*(f[i]+d[i]);
         }
         UB = 0.5*UB;
 */
         UB = UB - improv;
      }
                   
      /* mexPrintf("t=%d,k=%d, u=%d, tau1=%f, den1=%f, num1=%f, delta=%f\n", 
             t,k,u,tau1,den1,num1,delta);*/

    }

    /* -- Computing LB --------------------------------------*/

    /*
    LB = 0.5*x'*(f-d);   
    for k=1:n,
      LB = LB + b*min(d(find(I==k)));
    end */
    
    for( i=0, UB = 0, LB=0; i < n; i++) {
       UB += x[i]*(f[i]+d[i]);
       LB += x[i]*(f[i]-d[i]);
    }
    UB = 0.5*UB;
    LB = 0.5*LB;

    for( k=0; k < m; k++ ) { 
      for( j=0,tmp = PLUS_INF; j < nk[k]; j++ ) {
        i = inx[INDEX(j,k,n)];

/*      for( j=0, tmp = PLUS_INF; j < n; j++ ) {
        if( I[j]-1 == i ) tmp = MIN(tmp, d[j]);*/
        tmp = MIN(tmp, d[i]);
      }
      if( tmp < 0) LB += b*tmp;
    }

    /* Store LB and UB */
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

    if( verb > 0 && (exitflag > 0 || (t % verb)==0 )) {
       mexPrintf("%d: UB=%.10f, LB=%.10f, UB-LB=%.10f, (UB-LB)/|UB|=%.10f \n",
        t, UB, LB, UB-LB, (UB!=0) ? (UB-LB)/ABS(UB) : 0);      
    }    

  }

  /* -- Find which stopping consition has been used -------- */
  if( UB-LB < tolabs ) exitflag = 1;
  else if(UB-LB < ABS(UB)*tolrel ) exitflag = 2;
  else exitflag = 0;

  /*----------------------------------------------------------   
    Set up outputs                                          
  ---------------------------------------------------------- */
  (*ptr_t) = t;
  (*ptr_History) = History;

  /*----------------------------------------------------------
    Clean up
  ---------------------------------------------------------- */
  mxFree( d );
  mxFree( inx );
  mxFree( nk );
  
  return( exitflag ); 

}

