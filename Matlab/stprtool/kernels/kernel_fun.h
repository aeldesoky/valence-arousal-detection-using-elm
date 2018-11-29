/* --------------------------------------------------------------------
 Evaluation of kernel functions: linear, poly, rbf, sigmoidal.

 About: Statistical Pattern Recognition Toolbox
 (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
 <a href="http://www.cvut.cz">Czech Technical University Prague</a>
 <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
 <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

 Modifications:
 4-may-2004, VF
 14-November-2001, V.Franc, sigmoid kernel
 22-September-2001, V. Franc, created.
-------------------------------------------------------------------- */

/* --- Exported global variables ----------------------------------- */

extern double *dataA;   /* pointer at the first pattern */
extern double *dataB;   /* pointer at the second pattern */
extern long dim;        /* dimension of patterns */
extern int ker;         /* kernel type (0 - linear, 1 - polynom, 2 - rbf */
extern double *arg1;    /* argument of the kernel */
extern long ker_cnt;    /* number of kernel evaluations */
extern char *kernel_name[]; /* kernel names */


/* --- Exported functions ------------------------------------------ */

/* --------------------------------------------------------------------
 Computes kernel function for a-th (dataA) and b-th (dataB)
-------------------------------------------------------------------- */
extern double kernel( long a, long b );

/* --------------------------------------------------------------------
 Converts string kernel identifier to int.
-------------------------------------------------------------------- */
int kernel_id( const mxArray *prhs1 );
