/*-----------------------------------------------------------------------
gmnplib.c: Library of solvers for Generalized Minimal Norm Problem (GMNP).
 
 Generalized Minimal Norm Problem to solve is
  
  min 0.5*alpha'*H*alpha + c'*alpha

  subject to  sum(alpha) = 1,  alpha(i) >= 0
  
 H [dim x dim] is symmetric positive definite matrix.
 c [dim x 1] is an arbitrary vector.

 The precision of the found solution is given by
 the parameters tmax, tolabs and tolrel which
 define the stopping conditions:
 
 UB-LB <= tolabs      ->  exit_flag = 1   Abs. tolerance.
 UB-LB <= UB*tolrel   ->  exit_flag = 2   Relative tolerance.
 LB > th              ->  exit_flag = 3   Threshold on lower bound.
 t >= tmax            ->  exit_flag = 0   Number of iterations.

 UB ... Upper bound on the optimal solution.
 LB ... Lower bound on the optimal solution.
 t  ... Number of iterations.
 History ... Value of LB and UB wrt. number of iterations.


 The following algorithms are implemented:
 ..............................................

 - GMNP solver based on MDM algorithm.
   exitflag = gmnp_mdm( &get_col, diag_H, vector_c, dim,  
                 tmax, tolabs, tolrel, th, &alpha, &t, &History, verb  );

 - GMNP solver based on improved MDM algorithm 1 (u fixed v optimized)
    exitflag = gmnp_imdm( &get_col, diag_H, vector_c, dim,  
                 tmax, tolabs, tolrel, th, &alpha, &t, &History, verb  );

 - GMNP solver based on improved MDM algorithm 2 (u fixed v optimized 
     and vice versa)
    exitflag = gmnp_iimdm( &get_col, diag_H, vector_c, dim,  
                 tmax, tolabs, tolrel, th, &alpha, &t, &History, verb  );

 - GMNP solver based on the Kowalczyk's algorithm.
    exitflag = gmnp_kowalczyk( &get_col, diag_H, vector_c, dim,  
                  tmax, tolabs, tolrel, th, &alpha, &t, &History, verb  ); 

 - GMNP solver based on the Keerthis's algorithm.
    exitflag = gmnp_keerthi( &get_col, diag_H, vector_c, dim, 
                  tmax, tolabs, tolrel, th, &alpha, &t, &History, verb );

 - GMNP solver based on the Kozinec (alis Gilbert's) algorithm.
    exitflag = gmnp_kozinec( &get_col, diag_H, vector_c, dim, 
                  tmax, tolabs, tolrel, th, &alpha, &t, &History, verb );

  For more info refer to V.Franc: Optimization Algorithms for Kernel 
  Methods. Research report. CTU-CMP-2005-22. CTU FEL Prague. 2005.
  ftp://cmp.felk.cvut.cz/pub/cmp/articles/franc/Franc-PhD.pdf .

 Modifications:
 09-sep-2005, VF
 24-jan-2005, VF
 26-nov-2004, VF
 25-nov-2004, VF
 21-nov-2004, VF
 20-nov-2004, VF
 31-may-2004, VF
 23-Jan-2004, VF

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

#define ABS(A) ((A >= 0) ? A : -A)
#define MIN(A,B) ((A < B) ? A : B)
#define INDEX(ROW,COL,DIM) ((COL*DIM)+ROW)


/* --------------------------------------------------------------
 GMNP solver based on MDM algorithm.

 Usage: exitflag = gmnp_mdm( &get_col, diag_H, vector_c, dim,  
                  tmax, tolabs, tolrel, th, &alpha, &t, &History );
-------------------------------------------------------------- */
int gmnp_mdm(const void* (*get_col)(long,long),
            double *diag_H,
            double *vector_c,
            long dim, 
            long tmax,
            double tolabs,
            double tolrel,
            double th,
            double *alpha,
            long  *ptr_t,
            double **ptr_History,
            long verb)
{
  double LB;
  double UB;
  double aHa, ac;
  double tmp, tmp1;
  double Huu, Huv, Hvv;
  double min_beta, max_beta, beta;
  double lambda;
  double *History;
  double *Ha;
  double *tmp_ptr;
  double *col_u, *col_v;
  long u;
  long v;
  long new_u;
  long new_v;
  long i;
  long t;
  long History_size;
  int exitflag;

  /* ------------------------------------------------------------ */
  /* Initialization                                               */
  /* ------------------------------------------------------------ */

  Ha = mxCalloc(dim, sizeof(double));
  if( Ha == NULL ) mexErrMsgTxt("Not enough memory.");

  History_size = (tmax < HISTORY_BUF ) ? tmax+1 : HISTORY_BUF;
  History = mxCalloc(History_size*2,sizeof(double));
  if( History == NULL ) mexErrMsgTxt("Not enough memory.");

  /* inx = argmin(0.5*diag_H + vector_c ); */
  for( tmp1 =  PLUS_INF, i = 0; i < dim; i++ ) {
    tmp = 0.5*diag_H[i] + vector_c[i];
    if( tmp1 > tmp) {
      tmp1 = tmp;
      v = i;
    }
  }
  col_v = (double*)get_col(v,-1);

  for( min_beta = PLUS_INF, i = 0; i < dim; i++ ) 
  {
    alpha[i] = 0;
    Ha[i] = col_v[i];

    beta = Ha[i] + vector_c[i];
    if( beta < min_beta ) {
      min_beta = beta;
      u = i;
    }
  }

  alpha[v] = 1;
  aHa = diag_H[v];
  ac = vector_c[v];

  UB = 0.5*aHa + ac;
  LB = min_beta - 0.5*aHa;
  t = 0;
  History[INDEX(0,0,2)] = LB;
  History[INDEX(1,0,2)] = UB;

  if( verb ) {
    mexPrintf("Init: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
      UB, LB, UB-LB,(UB-LB)/UB);
  }  

  /* Stopping conditions */
  if( UB-LB <= tolabs ) exitflag = 1;
  else if(UB-LB <= ABS(UB)*tolrel ) exitflag = 2;
  else if(LB > th ) exitflag = 3;
  else exitflag = -1;

  /* ------------------------------------------------------------ */
  /* Main optimization loop                                       */
  /* ------------------------------------------------------------ */

  while( exitflag == -1 ) 
  {
    t++;     

    col_u = (double*)get_col(u,-1);
    col_v = (double*)get_col(v,u);

    /* Adaptation rule and update */
    Huu = diag_H[u];
    Hvv = diag_H[v];
    Huv = col_u[v];

    lambda = (Ha[v]-Ha[u]+vector_c[v]-vector_c[u])/(alpha[v]*(Huu-2*Huv+Hvv));
    if( lambda < 0 ) lambda = 0; else if (lambda > 1) lambda = 1;

    aHa = aHa + 2*alpha[v]*lambda*(Ha[u]-Ha[v])+
                lambda*lambda*alpha[v]*alpha[v]*(Huu-2*Huv+Hvv);

    ac = ac + lambda*alpha[v]*(vector_c[u]-vector_c[v]);

    tmp = alpha[v];
    alpha[u]=alpha[u]+lambda*alpha[v];
    alpha[v]=alpha[v]-lambda*alpha[v];

    UB = 0.5*aHa + ac;

    min_beta = PLUS_INF; 
    max_beta = MINUS_INF;
    for( i = 0; i < dim; i++ ) 
    {
       Ha[i] = Ha[i] + lambda*tmp*(col_u[i] - col_v[i]);

       beta = Ha[i]+ vector_c[i];

       if( alpha[i] !=0 && max_beta < beta ) 
       {
         new_v = i;
         max_beta = beta;
       }

       if( beta < min_beta )
       { 
         new_u = i;
         min_beta = beta;
       }
    }    

    LB = min_beta - 0.5*aHa; 
    u = new_u;    
    v = new_v;

    /* Stopping conditions */
    if( UB-LB <= tolabs ) exitflag = 1; 
    else if( UB-LB <= ABS(UB)*tolrel ) exitflag = 2;
    else if(LB > th ) exitflag = 3;
    else if(t >= tmax) exitflag = 0; 

    if( verb && (t % verb) == 0) {
      mexPrintf("%d: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
        t, UB, LB, UB-LB,(UB-LB)/UB);
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

  /* print info about last iteration*/
  if(verb && (t % verb) ) {
    mexPrintf("exit: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
      UB, LB, UB-LB,(UB-LB)/UB);
  }  


  /*------------------------------------------------------- */
  /* Set outputs                                            */
  /*------------------------------------------------------- */
  (*ptr_t) = t;
  (*ptr_History) = History;

  /* Free memory */
  mxFree( Ha );
  
  return( exitflag ); 
}


/* --------------------------------------------------------------
 GMNP solver based on improved MDM algorithm 1.

 Search strategy: u determined by common rule and v is 
 optimized.

 Usage: exitflag = gmnp_imdm( &get_col, diag_H, vector_c, dim,  
                  tmax, tolabs, tolrel, th, &alpha, &t, &History );
-------------------------------------------------------------- */

int gmnp_imdm(const void* (*get_col)(long,long),
            double *diag_H,
            double *vector_c,
            long dim, 
            long tmax,
            double tolabs,
            double tolrel,
            double th,
            double *alpha,
            long  *ptr_t,
            double **ptr_History,
            long verb)
{
  double LB;
  double UB;
  double aHa, ac;
  double tmp, tmp1;
  double Huu, Huv, Hvv;
  double min_beta, max_beta, beta;
  double max_improv, improv;
  double lambda;
  double *History;
  double *Ha;
  double *tmp_ptr;
  double *col_u, *col_v;
  long u;
  long v;
  long new_u;
  long new_v;
  long i;
  long t;
  long History_size;
  int exitflag;

  /* ------------------------------------------------------------ */
  /* Initialization                                               */
  /* ------------------------------------------------------------ */

  Ha = mxCalloc(dim, sizeof(double));
  if( Ha == NULL ) mexErrMsgTxt("Not enough memory.");

  History_size = (tmax < HISTORY_BUF ) ? tmax+1 : HISTORY_BUF;
  History = mxCalloc(History_size*2,sizeof(double));
  if( History == NULL ) mexErrMsgTxt("Not enough memory.");

  /* inx = argmin(0.5*diag_H + vector_c ); */
  for( tmp1 =  PLUS_INF, i = 0; i < dim; i++ ) {
    tmp = 0.5*diag_H[i] + vector_c[i];
    if( tmp1 > tmp) {
      tmp1 = tmp;
      v = i;
    }
  }

  col_v = (double*)get_col(v,-1);

  for( min_beta = PLUS_INF, i = 0; i < dim; i++ ) 
  {
    alpha[i] = 0;
    Ha[i] = col_v[i];

    beta = Ha[i] + vector_c[i];
    if( beta < min_beta ) {
      min_beta = beta;
      u = i;
    }
  }

  alpha[v] = 1;
  aHa = diag_H[v];
  ac = vector_c[v];

  UB = 0.5*aHa + ac;
  LB = min_beta - 0.5*aHa;

  t = 0;
  History[INDEX(0,0,2)] = LB;
  History[INDEX(1,0,2)] = UB;

  if( verb ) {
    mexPrintf("Init: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
      UB, LB, UB-LB,(UB-LB)/UB);
  }  

  /* Stopping conditions */
  if( UB-LB <= tolabs ) exitflag = 1;
  else if(UB-LB <= ABS(UB)*tolrel ) exitflag = 2;
  else if(LB > th ) exitflag = 3;
  else exitflag = -1;

  /* ------------------------------------------------------------ */
  /* Main optimization loop                                       */
  /* ------------------------------------------------------------ */

  col_u = (double*)get_col(u,-1);
  while( exitflag == -1 ) 
  {
    t++;     

    col_v = (double*)get_col(v,u);

    /* Adaptation rule and update */
    Huu = diag_H[u];
    Hvv = diag_H[v];
    Huv = col_u[v];

    lambda = (Ha[v]-Ha[u]+vector_c[v]-vector_c[u])/(alpha[v]*(Huu-2*Huv+Hvv));
    if( lambda < 0 ) lambda = 0; else if (lambda > 1) lambda = 1;

    aHa = aHa + 2*alpha[v]*lambda*(Ha[u]-Ha[v])+
                lambda*lambda*alpha[v]*alpha[v]*(Huu-2*Huv+Hvv);

    ac = ac + lambda*alpha[v]*(vector_c[u]-vector_c[v]);

    tmp = alpha[v];
    alpha[u]=alpha[u]+lambda*alpha[v];
    alpha[v]=alpha[v]-lambda*alpha[v];

    UB = 0.5*aHa + ac;
    
/*    max_beta = MINUS_INF;*/
    for( min_beta = PLUS_INF, i = 0; i < dim; i++ ) 
    {
       Ha[i] = Ha[i] + lambda*tmp*(col_u[i] - col_v[i]);

       beta = Ha[i]+ vector_c[i];

       if( beta < min_beta )
       { 
         new_u = i;
         min_beta = beta;
       }
    }    

    LB = min_beta - 0.5*aHa; 
    u = new_u;    
    col_u = (double*)get_col(u,-1);

    /* search for optimal v while u is fixed */
    for( max_improv =  MINUS_INF, i = 0; i < dim; i++ ) {

      if( alpha[i] != 0 ) {
        beta = Ha[i] + vector_c[i];

        if( beta >= min_beta ) {

          tmp = diag_H[u] - 2*col_u[i] + diag_H[i];
          if( tmp != 0 ) {
            improv = (0.5*(beta-min_beta)*(beta-min_beta))/tmp;

            if( improv > max_improv ) {
              max_improv = improv;
              v = i;
            }
          }
        }
      }
    }

    /* Stopping conditions */
    if( UB-LB <= tolabs ) exitflag = 1; 
    else if( UB-LB <= ABS(UB)*tolrel ) exitflag = 2;
    else if(LB > th ) exitflag = 3;
    else if(t >= tmax) exitflag = 0; 

    /* print info */
    if(verb && (t % verb) == 0 ) {
      mexPrintf("%d: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
        t, UB, LB, UB-LB,(UB-LB)/UB);
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

  /* print info about last iteration*/
  if(verb && (t % verb) ) {
    mexPrintf("exit: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
      UB, LB, UB-LB,(UB-LB)/UB);
  }  


  /*------------------------------------------------------- */
  /* Set outputs                                            */
  /*------------------------------------------------------- */
  (*ptr_t) = t;
  (*ptr_History) = History;

  /* Free memory */
  mxFree( Ha );
  
  return( exitflag ); 
}

/* --------------------------------------------------------------
 GMNP solver based on improved MDM algorithm 2.

 Search strategy: u fix and v optimzed plus v fixed and u 
 optimized. 

 Usage: exitflag = gmnp_iimdm( &get_col, diag_H, vector_c, dim,  
                  tmax, tolabs, tolrel, th, &alpha, &t, &History );
-------------------------------------------------------------- */

int gmnp_iimdm(const void* (*get_col)(long,long),
            double *diag_H,
            double *vector_c,
            long dim, 
            long tmax,
            double tolabs,
            double tolrel,
            double th,
            double *alpha,
            long  *ptr_t,
            double **ptr_History,
            long verb)
{
  double LB;
  double UB;
  double aHa, ac;
  double tmp, tmp1;
  double Huu, Huv, Hvv;
  double min_beta, max_beta, beta;
  double max_improv1, max_improv2, improv;
  double lambda;
  double *History;
  double *Ha;
  double *tmp_ptr;
  double *col_u, *col_v;
  long u;
  long v;
  long new_u;
  long new_v;
  long i;
  long t;
  long History_size;
  int exitflag;

  /* ------------------------------------------------------------ */
  /* Initialization                                               */
  /* ------------------------------------------------------------ */

  Ha = mxCalloc(dim, sizeof(double));
  if( Ha == NULL ) mexErrMsgTxt("Not enough memory.");

  History_size = (tmax < HISTORY_BUF ) ? tmax+1 : HISTORY_BUF;
  History = mxCalloc(History_size*2,sizeof(double));
  if( History == NULL ) mexErrMsgTxt("Not enough memory.");

  /* inx = argmin(0.5*diag_H + vector_c ); */
  for( tmp1 =  PLUS_INF, i = 0; i < dim; i++ ) {
    tmp = 0.5*diag_H[i] + vector_c[i];
    if( tmp1 > tmp) {
      tmp1 = tmp;
      v = i;
    }
  }

  col_v = (double*)get_col(v,-1);

  for( min_beta = PLUS_INF, i = 0; i < dim; i++ ) 
  {
    alpha[i] = 0;
    Ha[i] = col_v[i];

    beta = Ha[i] + vector_c[i];
    if( beta < min_beta ) {
      min_beta = beta;
      u = i;
    }
  }

  alpha[v] = 1;
  aHa = diag_H[v];
  ac = vector_c[v];

  UB = 0.5*aHa + ac;
  LB = min_beta - 0.5*aHa;

  t = 0;
  History[INDEX(0,0,2)] = LB;
  History[INDEX(1,0,2)] = UB;

  if( verb ) {
    mexPrintf("Init: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
      UB, LB, UB-LB,(UB-LB)/UB);
  }  

  /* Stopping conditions */
  if( UB-LB <= tolabs ) exitflag = 1;
  else if(UB-LB <= ABS(UB)*tolrel ) exitflag = 2;
  else if(LB > th ) exitflag = 3;
  else exitflag = -1;

  /* ------------------------------------------------------------ */
  /* Main optimization loop                                       */
  /* ------------------------------------------------------------ */

  col_u = (double*)get_col(u,-1);
  col_v = (double*)get_col(v,u);
  while( exitflag == -1 ) 
  {
    t++;     

    /* Adaptation rule and update */
    Huu = diag_H[u];
    Hvv = diag_H[v];
    Huv = col_u[v];

    lambda = (Ha[v]-Ha[u]+vector_c[v]-vector_c[u])/(alpha[v]*(Huu-2*Huv+Hvv));
    if( lambda < 0 ) lambda = 0; else if (lambda > 1) lambda = 1;

    aHa = aHa + 2*alpha[v]*lambda*(Ha[u]-Ha[v])+
                lambda*lambda*alpha[v]*alpha[v]*(Huu-2*Huv+Hvv);

    ac = ac + lambda*alpha[v]*(vector_c[u]-vector_c[v]);

    tmp = alpha[v];
    alpha[u]=alpha[u]+lambda*alpha[v];
    alpha[v]=alpha[v]-lambda*alpha[v];

    UB = 0.5*aHa + ac;
    
    for( max_beta = MINUS_INF, min_beta = PLUS_INF, i = 0; i < dim; i++ ) 
    {
       Ha[i] = Ha[i] + lambda*tmp*(col_u[i] - col_v[i]);

       beta = Ha[i]+ vector_c[i];

       if( beta < min_beta )
       { 
         new_u = i;
         min_beta = beta;
       }

       if( alpha[i] != 0 && max_beta < beta ) 
       {
         new_v = i;
         max_beta = beta;
       }
    }    

    LB = min_beta - 0.5*aHa; 

    col_u = (double*)get_col(new_u,-1);
    col_v = (double*)get_col(new_v,new_u);

    /* search for optimal v while u is fixed */
    max_improv1 =  MINUS_INF; max_improv2 =  MINUS_INF;
    for( i = 0; i < dim; i++ ) {

      beta = Ha[i] + vector_c[i]; 

      if( alpha[i] != 0 && beta > min_beta ) {

        tmp = diag_H[new_u] - 2*col_u[i] + diag_H[i];
        if( tmp != 0 ) 
        { 

          if((beta-min_beta)/(alpha[i]*tmp) < 1) 
          {
            improv = (0.5*(beta-min_beta)*(beta-min_beta))/tmp; 
          } else {
            improv = alpha[i]*(beta-min_beta) - 0.5*alpha[i]*alpha[i]*tmp;
          }

          if( improv > max_improv1 ) 
          {
            max_improv1 = improv;
            v = i;
          }
        }
      }

      if( max_beta > beta ) {
        tmp = diag_H[new_v] - 2*col_v[i] + diag_H[i];
        if( tmp != 0 ) 
        { 
          
          if((max_beta-beta)/(alpha[new_v]*tmp) < 1 ) 
          {
            improv = (0.5*(max_beta-beta)*(max_beta-beta))/tmp; 
          } else {
            improv = alpha[new_v]*(max_beta-beta)
                 -0.5*alpha[new_v]*alpha[new_v]*tmp;
          }

          if( improv > max_improv2 ) 
          {
            max_improv2 = improv;
            u = i;
          }
        }
      }
    }

    if( max_improv1 > max_improv2 ) {
      u = new_u;
      col_v = (double*)get_col(v,u);
    } else {
      v = new_v;
      col_u = (double*)get_col(u,v);
    }

/*    tmp = Ha[v] + vector_c[v] - Ha[u] - vector_c[u];*/

    /* Stopping conditions */
    if( UB-LB <= tolabs ) exitflag = 1; 
    else if( UB-LB <= ABS(UB)*tolrel ) exitflag = 2;
    else if(LB > th ) exitflag = 3;
    else if(t >= tmax) exitflag = 0; 

    /* print info */
    if(verb && (t % verb) == 0 ) {
      mexPrintf("%d: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
        t, UB, LB, UB-LB,(UB-LB)/UB);
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

  /* print info about last iteration*/
  if(verb && (t % verb) ) {
    mexPrintf("exit: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
      UB, LB, UB-LB,(UB-LB)/UB);
  }  

  /*------------------------------------------------------- */
  /* Set outputs                                            */
  /*------------------------------------------------------- */
  (*ptr_t) = t;
  (*ptr_History) = History;

  /* Free memory */
  mxFree( Ha );
  
  return( exitflag ); 
}


/* --------------------------------------------------------------
 GMNP solver based on Keerthi's algorithm.

 Usage: exitflag = gmnp_keerthi( &get_col, diag_H, vector_c, dim, 
                  tmax, tolabs, tolrel, th, &alpha, &t, &History );
-------------------------------------------------------------- */
int gmnp_keerthi(const void* (*get_col)(long,long),
            double *diag_H,
            double *vector_c,
            long dim, 
            long tmax,
            double tolabs,
            double tolrel,
            double th,
            double *alpha,
            long  *ptr_t,
            double **ptr_History,
            long verb)
{
  double LB;
  double UB;
  double aHa, ac;
  double den, tmp, tmp1;
  double Huu, Huv, Hvv;
  double min_beta, max_beta, beta;
  double gamma, omega;
  double *History;
  double *Ha;
  double *tmp_ptr;
  double x11, x12, x13, x23, x22, x33, x10, x20, x30;
  double a1, a2, a3, a4, a5;
  double UB123, UB1, UB2, UB3;
  double gamma1, gamma2, gamma3;
  double tmp_aHa1, tmp_aHa2, tmp_aHa3;
  double tmp_ac1, tmp_ac2, tmp_ac3;
  double *col_u, *col_v;
  long u;
  long v;
  long i;
  long t;
  long History_size;
  int exitflag;
  int nearest_segment;

  /* ------------------------------------------------------------ */
  /* Initialization                                               */
  /* ------------------------------------------------------------ */

  Ha = mxCalloc(dim, sizeof(double));
  if( Ha == NULL ) mexErrMsgTxt("Not enough memory.");

  History_size = (tmax < HISTORY_BUF ) ? tmax+1 : HISTORY_BUF;
  History = mxCalloc(History_size*2,sizeof(double));
  if( History == NULL ) mexErrMsgTxt("Not enough memory.");

  /* inx = argmin(0.5*diag_H + vector_c ); */
  for( tmp1 =  PLUS_INF, i = 0; i < dim; i++ ) {
    tmp = 0.5*diag_H[i] + vector_c[i];
    if( tmp1 > tmp) {
      tmp1 = tmp;
      v = i;
    }
  }

  col_v = (double*)get_col(v,-1);

  for( min_beta = PLUS_INF, i = 0; i < dim; i++ ) 
  {
    alpha[i] = 0;
    Ha[i] = col_v[i];

    beta = Ha[i] + vector_c[i];
    if( beta < min_beta ) {
      min_beta = beta;
      u = i;
    }
  }

  alpha[v] = 1;
  aHa = diag_H[v];
  ac = vector_c[v];

  UB = 0.5*aHa + ac;
  LB = min_beta - 0.5*aHa;
  t = 0;
  History[INDEX(0,0,2)] = LB;
  History[INDEX(1,0,2)] = UB;

  if( verb ) {
    mexPrintf("Init: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
      UB, LB, UB-LB,(UB-LB)/UB);
  }  

  /* Stopping conditions */
  if( UB-LB <= tolabs ) exitflag = 1;
  else if( UB-LB <= ABS(UB)*tolrel ) exitflag = 2;
  else if( LB > th ) exitflag = 3;
  else exitflag = -1;

  /* ------------------------------------------------------------ */
  /* Main optimization loop                                       */
  /* ------------------------------------------------------------ */

  while( exitflag == -1 ) 
  {
    t++;     

    col_u = (double*)get_col(u,-1);
    col_v = (double*)get_col(v,u);

    /* Adaptation rule and update */
    Huu = diag_H[u];
    Hvv = diag_H[v];
    Huv = col_u[v];

    x11 = aHa;
    x12 = Ha[u];
    x13 = aHa + alpha[v]*(Ha[u]-Ha[v]);
    x22 = Huu;
    x23 = Ha[u] + alpha[v]*(Huu-Huv);
    x33 = aHa + 2*alpha[v]*(Ha[u]-Ha[v])+alpha[v]*alpha[v]*(Huu-2*Huv+Hvv);

    x10 = ac;
    x20 = vector_c[u];
    x30 = ac + alpha[v]*(vector_c[u]-vector_c[v]);

    a1 = x11 - x12 - x13 + x23;
    a2 = x11 - 2*x12 + x22;
    a3 = x12 - x11 + x20 - x10;
    a4 = x11 - 2*x13 + x33;
    a5 = x13 - x11 + x30 - x10;

    den = a1*a1 - a2*a4;
    if( den ) {

      gamma = (a3*a4-a1*a5)/den;
      omega = (a2*a5-a3*a1)/den;

      if( gamma > 0 && omega > 0 && 1-gamma-omega > 0 ) {

        /* Ha = Ha*(1-gamma) + H(:,u)*(gamma+alpha(v)*omega)-H(:,v)*alpha(v)*omega;*/
        tmp = alpha[v]*omega;
        for(i = 0; i < dim; i++ ) {
          Ha[i] = Ha[i]*(1-gamma) + col_u[i]*(gamma+tmp) - col_v[i]*tmp;
        }

        /* 
          aHa = (1-omega-gamma)^2*x11 + gamma^2*x22 + omega^2*x33 + ...
          2*(1-omega-gamma)*gamma*x12 + 2*(1-omega-gamma)*omega*x13 + ...
          2*gamma*omega*x23; 
        */
        aHa = (1-omega-gamma)*(1-omega-gamma)*x11 + gamma*gamma*x22
               + omega*omega*x33 + 2*(1-omega-gamma)*gamma*x12 
               + 2*(1-omega-gamma)*omega*x13 + 2*gamma*omega*x23;
        ac = (1-gamma-omega)*x10 + gamma*x20 + omega*x30;


        /*
          alpha1 = zeros(dim,1);
          alpha1(u) = 1;
          alpha2 = alpha;
          alpha2(u) = alpha(u)+alpha(v);
          alpha2(v) = 0;
          alpha = alpha*(1-gamma-omega) + alpha1*gamma + alpha2*omega;
        */
        for(i = 0; i < dim; i++ ) {
          alpha[i] = alpha[i]*(1-gamma);
        }
       
        alpha[u] = alpha[u] + gamma + tmp;
        alpha[v] = alpha[v] - tmp;
          
        UB123 = 0.5*aHa + ac;
      } 
      else
      {
        UB123 = PLUS_INF;
      }

    }
    else
    {
      UB123 = PLUS_INF;
    }

    if( UB123 == PLUS_INF ) {
      /*
        % line segment between alpha and alpha1
        gamma1 = (x11-x12+x10-x20)/(x11-2*x12+x22);
        gamma1 = min(1,gamma1);
        tmp_aHa1 = (1-gamma1)^2*x11+2*gamma1*(1-gamma1)*x12+gamma1^2*x22;
        tmp_ac1 = (1-gamma1)*x10 + gamma1*x20; 
        UB1 = 0.5*tmp_aHa1 + tmp_ac1;
      */
      gamma1 = (x11-x12+x10-x20)/(x11-2*x12+x22);
      gamma1 = MIN(1,gamma1);
      tmp_aHa1 = (1-gamma1)*(1-gamma1)*x11+2*gamma1*(1-gamma1)*x12+gamma1*gamma1*x22;
      tmp_ac1 = (1-gamma1)*x10 + gamma1*x20; 
      UB1 = 0.5*tmp_aHa1 + tmp_ac1;
      
      /*
        % line segment between alpha and alpha2
        gamma2 = (x11-x13+x10-x30)/(x11-2*x13+x33);
        gamma2 = min(1,gamma2);
        tmp_aHa2 = (1-gamma2)^2*x11+2*gamma2*(1-gamma2)*x13+gamma2^2*x33;
        tmp_ac2 = (1-gamma2)*x10 + gamma2*x30;
        UB2 = 0.5*tmp_aHa2 + tmp_ac2;
      */
      gamma2 = (x11-x13+x10-x30)/(x11-2*x13+x33);
      gamma2 = MIN(1,gamma2);
      tmp_aHa2 = (1-gamma2)*(1-gamma2)*x11+2*gamma2*(1-gamma2)*x13+gamma2*gamma2*x33;
      tmp_ac2 = (1-gamma2)*x10 + gamma2*x30;
      UB2 = 0.5*tmp_aHa2 + tmp_ac2;
 
      /*
        % line segment between alpha1 and alpha2
        tmp_den = (x22-2*x23+x33);
        if tmp_den ~= 0,
          gamma3 = (x22-x23+x20-x30)/tmp_den;
          if gamma3 > 1 gamma3 = 1; end
          if gamma3 < 0 gamma3 = 0; end
          tmp_aHa3 = (1-gamma3)^2*x22+2*gamma3*(1-gamma3)*x23+gamma3^2*x33;
          tmp_ac3 = (1-gamma3)*x20 + gamma3*x30;
          UB3 = 0.5*tmp_aHa3 + tmp_ac3;
        else
          UB3 = UB;
        end
      */
      den = (x22-2*x23+x33);
      if( den ) {
        gamma3 = (x22-x23+x20-x30)/den;
        if( gamma3 > 1) gamma3 = 1; 
        if( gamma3 < 0) gamma3 = 0; 
        tmp_aHa3 = (1-gamma3)*(1-gamma3)*x22+2*gamma3*(1-gamma3)*x23+
           gamma3*gamma3*x33;
        tmp_ac3 = (1-gamma3)*x20 + gamma3*x30;
        UB3 = 0.5*tmp_aHa3 + tmp_ac3;
      }
      else
      {
        UB3 = UB;
      }

      /* nearest_segment = argmin( UB1, UB2, UB3 ) */
      if( UB1 <= UB2 ) {
        if( UB1 <= UB3 ) nearest_segment = 1; else nearest_segment = 3; 
      } 
      else 
      {
        if( UB2 <= UB3 ) nearest_segment = 2; else nearest_segment = 3;
      }


        /*
          alpha1 = zeros(dim,1);
          alpha1(u) = 1;
          alpha2 = alpha;
          alpha2(u) = alpha(u)+alpha(v);
          alpha2(v) = 0;
        */

      switch( nearest_segment )
      {
        case 1: 
          aHa = tmp_aHa1;
          ac = tmp_ac1;
          /*  
            Ha = Ha*(1-gamma1) + gamma1*H(:,u); 
            alpha = alpha*(1-gamma1)+gamma1*alpha1;
          */
          for( i = 0; i < dim; i++ ) {
            Ha[i] = Ha[i]*(1-gamma1) + gamma1*col_u[i];
            alpha[i] = alpha[i]*(1-gamma1);
          }
          alpha[u] += gamma1;
          break;

        case 2: 
          aHa = tmp_aHa2;
          ac = tmp_ac2;
          /*
            Ha = Ha + gamma2*alpha(v)*(H(:,u)-H(:,v));
            alpha = alpha*(1-gamma2)+gamma2*alpha2;
          */
          tmp = alpha[v]*gamma2;
          for( i = 0; i < dim; i++ ) {
            Ha[i] = Ha[i] + tmp*(col_u[i] - col_v[i]);
          }
          alpha[u] = alpha[u] + tmp;
          alpha[v] = alpha[v] - tmp;
          break;

        case 3: 
          aHa = tmp_aHa3;
          ac = tmp_ac3;
          /* 
            Ha = gamma3*Ha + H(:,u)*(1-gamma3+gamma3*alpha(v)) - ...
                 gamma3*alpha(v)*H(:,v);
            alpha = alpha1*(1-gamma3)+gamma3*alpha2;
          */
          tmp = alpha[v]*gamma3;
          for( i = 0; i < dim; i++ ) {
            Ha[i] = gamma3*Ha[i] + col_u[i]*(1-gamma3+tmp) - tmp*col_v[i];
            alpha[i] = alpha[i]*gamma3;
          }          
          alpha[u] = alpha[u] + 1 - gamma3 + tmp;
          alpha[v] = alpha[v] - tmp;

          break;
      }  
    }

    UB = 0.5*aHa + ac;

    min_beta = PLUS_INF; 
    max_beta = MINUS_INF;
    for( i = 0; i < dim; i++ ) 
    {
       beta = Ha[i]+ vector_c[i];

       if( alpha[i] !=0 && max_beta < beta ) 
       {
         v = i;
         max_beta = beta;
       }

       if( beta < min_beta )
       { 
         u = i;
         min_beta = beta;
       }
    }    

    LB = min_beta - 0.5*aHa; 

    /* Stopping conditions */
    if( UB-LB <= tolabs ) exitflag = 1; 
    else if( UB-LB <= ABS(UB)*tolrel ) exitflag = 2;
    else if(LB > th ) exitflag = 3;
    else if(t >= tmax) exitflag = 0; 

    /* print info */
    if(verb && (t % verb) == 0 ) {
      mexPrintf("%d: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
        t, UB, LB, UB-LB,(UB-LB)/UB);
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

  /* print info about last iteration*/
  if(verb && (t % verb) ) {
    mexPrintf("exit: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
      UB, LB, UB-LB,(UB-LB)/UB);
  }  

  /*------------------------------------------------------- */
  /* Set outputs                                            */
  /*------------------------------------------------------- */
  (*ptr_t) = t;
  (*ptr_History) = History;

  /* Free memory */
  mxFree( Ha );
  
  return( exitflag ); 
}


/* --------------------------------------------------------------
 GMNP solver based on the Kowalczyk's algorithm 
 (Maximal Margin Perceptron).

 Usage: exitflag = gmnp_kowalczyk( &get_col, diag_H, vector_c, dim,  
                  tmax, tolabs, tolrel, th, &alpha, &t, &History ); 
-------------------------------------------------------------- */
int gmnp_kowalczyk(const void* (*get_col)(long,long),
            double *diag_H,
            double *vector_c,
            long dim, 
            long tmax,
            double tolabs,
            double tolrel,
            double th,
            double *alpha,
            long  *ptr_t,
            double **ptr_History,
            long verb)
{
  double LB;
  double UB;
  double aHa, ac;
  double tmp, tmp1;
  double tmp_aHa, tmp_ac, tmp_UB, tmp_gamma;
  double min_beta, beta;
  double x10, x20, x11, x12, x22;
  double delta;
  double gamma;
  double *History;
  double *Ha;
  double *tmp_ptr;
  double *col_inx;
  long i;
  long inx;
  long t;
  long History_size;
  int exitflag;

  /* ------------------------------------------------------------ */
  /* Initialization                                               */
  /* ------------------------------------------------------------ */

  Ha = mxCalloc(dim, sizeof(double));
  if( Ha == NULL ) mexErrMsgTxt("Not enough memory.");

  History_size = (tmax < HISTORY_BUF ) ? tmax+1 : HISTORY_BUF;
  History = mxCalloc(History_size*2,sizeof(double));
  if( History == NULL ) mexErrMsgTxt("Not enough memory.");

  /* inx = argmin(0.5*diag_H + vector_c ); */
  for( tmp1 =  PLUS_INF, i = 0; i < dim; i++ ) {
    tmp = 0.5*diag_H[i] + vector_c[i];
    if( tmp1 > tmp) {
      tmp1 = tmp;
      inx = i;
    }
  }

  col_inx = (double*)get_col(inx,-1);

  for( min_beta = PLUS_INF, i = 0; i < dim; i++ ) 
  {
    alpha[i] = 0;
    Ha[i] = col_inx[i];

    beta = Ha[i] + vector_c[i];
    if( beta < min_beta ) {
      min_beta = beta;
    }
  }

  alpha[inx] = 1;
  aHa = diag_H[inx];
  ac = vector_c[inx];

  UB = 0.5*aHa + ac;
  LB = min_beta - 0.5*aHa;
  t = 0;
  History[INDEX(0,0,2)] = LB;
  History[INDEX(1,0,2)] = UB;

  if( verb ) {
    mexPrintf("Init: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
      UB, LB, UB-LB,(UB-LB)/UB);
  }  

  /* Stopping conditions */
  if( UB-LB <= tolabs ) exitflag = 1;
  else if(UB-LB <= ABS(UB)*tolrel ) exitflag = 2;
  else if(LB > th ) exitflag = 3;
  else exitflag = -1;

  /* ------------------------------------------------------------ */
  /* Main optimization loop                                       */
  /* ------------------------------------------------------------ */

  while( exitflag == -1 ) 
  {
    t++;     

    x11 = aHa;
    x10 = ac;
    /* searches for rule which yileds the biggest improvement */
    for( i = 0; i < dim; i++ ) 
    {
       delta = Ha[i] + vector_c[i] - aHa - ac;
    
       tmp_UB = PLUS_INF;
       if( delta < 0 ) 
       {
          /* Kozinec rule */
          x12 = Ha[i];
          x20 = vector_c[i];
          x22 = diag_H[i];
      
          tmp_gamma = (x11 - x12 + x10 - x20)/(x11 - 2*x12 + x22);
          tmp_gamma = MIN(1,tmp_gamma);

          tmp_aHa = (1-tmp_gamma)*(1-tmp_gamma)*x11+2*(1-tmp_gamma)*tmp_gamma*x12
                  +tmp_gamma*tmp_gamma*x22;
          tmp_ac = (1-tmp_gamma)*x10+tmp_gamma*x20;
          tmp_UB = 0.5*tmp_aHa + tmp_ac;
       }
       else if( delta > 0 && alpha[i] < 1 && alpha[i] > 0)
       {
          x12 = (x11-alpha[i]*Ha[i])/(1-alpha[i]);
          x22 = (x11-2*alpha[i]*Ha[i]+alpha[i]*alpha[i]*diag_H[i])/
                 ((1-alpha[i])*(1-alpha[i]));
          x20 = (x10-alpha[i]*vector_c[i])/(1-alpha[i]);

          tmp_gamma = (x11 - x12 + x10 - x20)/(x11 - 2*x12 + x22);
          tmp_gamma = MIN(1,tmp_gamma);
  
          tmp_aHa = (1-tmp_gamma)*(1-tmp_gamma)*x11+2*(1-tmp_gamma)*tmp_gamma*x12
                  + tmp_gamma*tmp_gamma*x22;
          tmp_ac = (1-tmp_gamma)*x10+tmp_gamma*x20;
          tmp_UB = 0.5*tmp_aHa + tmp_ac;
       }

       if( tmp_UB < UB ) 
       {
         UB = tmp_UB;
         gamma = tmp_gamma;
         aHa = tmp_aHa;
         ac = tmp_ac;
         inx = i;
       }
    }

    col_inx = (double*)get_col(inx,-1);

    /* Use the update with biggest improvement */
    delta = Ha[inx] + vector_c[inx] - x11 - x10;
    if( delta < 0 ) 
    {
       /* Kozinec rule */
      for(i = 0; i < dim; i++ ) {
        Ha[i] = Ha[i]*(1-gamma) + gamma*col_inx[i];
        alpha[i] = alpha[i]*(1-gamma);
      }
      alpha[inx] = alpha[inx] + gamma;
    }
    else
    {
      /* Inverse Kozinec rule */
      tmp = gamma*alpha[inx];
      tmp1 = 1-alpha[inx];
      for(i = 0; i < dim; i++ ) {
        Ha[i] = (Ha[i]*(tmp1+tmp) - tmp*col_inx[i])/tmp1; 

        alpha[i] = alpha[i]*(1-gamma) + gamma*alpha[i]/tmp1;
      }
      alpha[inx] = alpha[inx] - tmp/tmp1;
    }

    min_beta = PLUS_INF; 
    for( i = 0; i < dim; i++ ) 
    {
       beta = Ha[i]+ vector_c[i];
       if( beta < min_beta ) { min_beta = beta; }
    }    

    LB = min_beta - 0.5*aHa; 

    /* Stopping conditions */
    if( UB-LB <= tolabs ) exitflag = 1; 
    else if( UB-LB <= ABS(UB)*tolrel ) exitflag = 2;
    else if( LB > th ) exitflag = 3;
    else if(t >= tmax) exitflag = 0; 

    /* print info */
    if(verb && (t % verb) == 0 ) {
      mexPrintf("%d: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
        t, UB, LB, UB-LB,(UB-LB)/UB);
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

  /* print info about last iteration*/
  if(verb && (t % verb) ) {
    mexPrintf("exit: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
      UB, LB, UB-LB,(UB-LB)/UB);
  }  

  /*------------------------------------------------------- */
  /* Set outputs                                            */
  /*------------------------------------------------------- */
  (*ptr_t) = t;
  (*ptr_History) = History;

  /* Free memory */
  mxFree( Ha );
  
  return( exitflag ); 
}

/* --------------------------------------------------------------
 GMNP solver based on the Kozinec algorithm.

 Usage: exitflag = gmnp_kozinec( &get_col, diag_H, vector_c, dim,  
                  tmax, tolabs, tolrel, th, &alpha, &t, &History );
-------------------------------------------------------------- */
int gmnp_kozinec(const void* (*get_col)(long,long),
            double *diag_H,
            double *vector_c,
            long dim, 
            long tmax,
            double tolabs,
            double tolrel,
            double th,
            double *alpha,
            long  *ptr_t,
            double **ptr_History,
            long verb)
{
  double LB;
  double UB;
  double aHa, ac;
  double tmp, tmp1;
  double min_beta, beta;
  double lambda;
  double *History;
  double *Ha;
  double *tmp_ptr;
  double *col_u;
  long u, inx;
  long new_u;
  long i;
  long t;
  long History_size;
  int exitflag;

  /* ------------------------------------------------------------ */
  /* Initialization                                               */
  /* ------------------------------------------------------------ */

  Ha = mxCalloc(dim, sizeof(double));
  if( Ha == NULL ) mexErrMsgTxt("Not enough memory.");

  History_size = (tmax < HISTORY_BUF ) ? tmax+1 : HISTORY_BUF;
  History = mxCalloc(History_size*2,sizeof(double));
  if( History == NULL ) mexErrMsgTxt("Not enough memory.");

  /* inx = argmin(0.5*diag_H + vector_c ); */
  for( tmp1 =  PLUS_INF, i = 0; i < dim; i++ ) {
    tmp = 0.5*diag_H[i] + vector_c[i];
    if( tmp1 > tmp) {
      tmp1 = tmp;
      inx = i;
    }
  }
  col_u = (double*)get_col(inx,-1);

  for( min_beta = PLUS_INF, i = 0; i < dim; i++ ) 
  {
    alpha[i] = 0;
    Ha[i] = col_u[i];

    beta = Ha[i] + vector_c[i];
    if( beta < min_beta ) {
      min_beta = beta;
      u = i;
    }
  }

  alpha[inx] = 1;
  aHa = diag_H[inx];
  ac = vector_c[inx];

  UB = 0.5*aHa + ac;
  LB = min_beta - 0.5*aHa;
  t = 0;
  History[INDEX(0,0,2)] = LB;
  History[INDEX(1,0,2)] = UB;

  if( verb ) {
    mexPrintf("Init: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
      UB, LB, UB-LB,(UB-LB)/UB);
  }  

  /* Stopping conditions */
  if( UB-LB <= tolabs ) exitflag = 1;
  else if(UB-LB <= ABS(UB)*tolrel ) exitflag = 2;
  else if(LB > th ) exitflag = 3;
  else exitflag = -1;

  /* ------------------------------------------------------------ */
  /* Main optimization loop                                       */
  /* ------------------------------------------------------------ */

  while( exitflag == -1 ) 
  {
    t++;     

    col_u = (double*)get_col(u,-1);

    /* Adaptation rule and update */
    lambda = (aHa - Ha[u] + ac - vector_c[u])/(aHa - 2*Ha[u] + diag_H[u]);
    lambda = MIN(1,lambda);

    aHa = aHa*(1-lambda)*(1-lambda) + 2*lambda*(1-lambda)*Ha[u] 
         + diag_H[u]*lambda*lambda;
    ac = ac*(1-lambda) + lambda*vector_c[u];

    min_beta = PLUS_INF;
    for( i = 0; i < dim; i++ ) {
      alpha[i] = alpha[i]*(1-lambda);
      Ha[i] = Ha[i]*(1-lambda) + lambda*col_u[i];

      beta = Ha[i] + vector_c[i];
      if( min_beta > beta ) { 
        min_beta = beta; 
        new_u = i; 
      }
    }
    alpha[u] = alpha[u] + lambda;

    UB = 0.5*aHa + ac;
    LB = min_beta - 0.5*aHa; 
    u = new_u;    


    /* Stopping conditions */
    if( UB-LB <= tolabs ) exitflag = 1; 
    else if( UB-LB <= ABS(UB)*tolrel ) exitflag = 2;
    else if( LB > th ) exitflag = 3;
    else if(t >= tmax) exitflag = 0; 

    /* print info about last iteration*/
    if(verb && (t % verb) == 0 ) {
      mexPrintf("%d: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
        t, UB, LB, UB-LB,(UB-LB)/UB);
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

  /* print info about last iteration*/
  if(verb && (t % verb) ) {
    mexPrintf("exit: UB=%f, LB=%f, UB-LB=%f, (UB-LB)/|UB|=%f \n",
      UB, LB, UB-LB,(UB-LB)/UB);
  }  

  /*------------------------------------------------------- */
  /* Set outputs                                            */
  /*------------------------------------------------------- */
  (*ptr_t) = t;
  (*ptr_History) = History;

  /* Free memory */
  mxFree( Ha );
  
  return( exitflag ); 
}

