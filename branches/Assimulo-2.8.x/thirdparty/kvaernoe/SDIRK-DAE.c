#include <stdio.h> 
#include <math.h> 
#include <stdlib.h> 
#include <float.h>

/* Definition of macros and constants. */ 
/*-------------------------------------*/ 
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y)) 
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

#define HMIN 1.e-8
#define EPS sqrt(DBL_EPSILON)
#define ROCK 0.5          /* Max allowed rate_of_convergence.           */
#define PES_FAC 0.9       /* Pessismist factor, used for stepsize       */
#define KAPPA 0.5         /* Iteration error < KAPPA*tolerance.         */

/* problem :                                                           */
/*---------------------------------------------------------------------*/ 
struct problem 
{  int n;               /* The dimension of the system.                */
   double *v;           /* The dependent variables.                    */
   double *vp;          /* The derivatives of the dependent variables. */ 
   double x;            /* The independent variable.                   */
};

/* information :                                                       */
/*---------------------------------------------------------------------*/ 
struct information 
{  int num_jacobian;               /* 0 if exact jacobian is supplied, */ 
                                   /* 1 if not.                        */ 
   double *abserr;                 /* Absolute and                     */ 
   double relerr;                  /* relative tolerance.              */ 
   double xend;                    /* End point for the integration.   */ 
   double h0;                      /* Starting stepsize.               */ 
   double hmax;                    /* Max allowed stepsize.            */
};

/* etcetera: work arrays and variables.                                    */
/*-------------------------------------------------------------------------*/
struct etcetera
{
   double *dfdv;                    /* The Jacobi matrix df/dv.            */
   double *dfdvp;                   /* The Jacobi matrix df/dv'.           */
   double *newton;                  /* The (factorised) Newton matrix.     */
   int *pivot;                      /* The pivot elements to newton.       */
   double *tolerance;               /* The array : abserr(i)+relerr*v(i).  */ 
   double min_tolerance;            /* Used for the residual error test.   */ 
   double *b;                       /* The b in f(Vi',h*g*VP+b,xi) = O.    */
   double rate_of_convergence;
   double h;                        /* The stepsize.                       */ 
   double gh;                       /* gh = a(i,i)*stepsize.               */ 
   double *temp1;                   /* Work arrays of dimension n.         */
   double *temp2;
   int jacobi_status;               /* 0 if the jacobian is new, else 1.   */ 
   int newton_status;               /* 0 if the Newton matrix is fact.     */
                                    /* with current stepsize, else 1.      */
   int new_jacobian;                /* Updated jacobians required.         */
   int new_fact;                    /* New LU-factorisation is required.   */ 
   int steps_red;                   /* Number of steps before the          */ 
                                    /* stepsize is allowed to increase.    */ 
   double *vp_stage[4];             /* Stage values. */
};


/* method : The parameters of the method.                           */
/*------------------------------------------------------------------*/
struct method
{
   int s;             /* Number of stages.                          */
   double a[4][4];
   double b_error[4]; /* The b-vector for error estimate.           */
   double c[4];
};

/* statistic :                                                      */
/*------------------------------------------------------------------*/
struct statistic
{
   int steps;
   int rejected_steps;
   int notkonv_steps;
   int function_evaluations; 
   int jacobi_evaluations;
   int factorisations;         /* of the Newton matrix.             */
   int solutions;              /* of the linear system.             */
} number_of;
/************ End of the header. ************************************/



/*------------------------------------------------------------------*/
/* sdirk_dae : Solves fully implicit uniform index 1 DAEs with a    */
/*             SDIRK method of order 3, embedded with a second      */
/*             order error estimating method.                       */
/*------------------------------------------------------------------*/

int sdirk_dae(struct problem *prob, struct information *info)

/*------------------------------------------------------------------*/
/* var : On entry : Function values at the starting point.          */
/*        On exit : Function values at the end point.               */
/* info : User supplied information. Not altered by the function.   */
/*                                                                  */
/* returns : 0 :The integration has been successfully performed.    */
/*           1 : The Newton matrix is singular.                     */
/*           2 : The Newton-iterations fails to converge.           */
/*           3 : The local error criteria fails to be satisfied,    */
/*               the step-size becomes too small.                   */
/*           4 : The error tolerance is too stringent.              */
/*------------------------------------------------------------------*/
{
   int i, flag;
   double step_fact;
   struct method rk;
   struct problem prob_old, *var_old, *var;
   struct etcetera etc;
   extern struct statistic number_of;

   /* Allocation of memory, and initialization of variables. */ 
   /*--------------------------------------------------------*/
   work_alloc(prob->n, &prob_old, &etc);
   var_old = &prob_old;
   var = prob;
   var_old->n = var->n;
   for (i = 0; i < var->n; ++i)
   {
      var_old->v[i] = var->v[i]; 
      var_old->vp[i] = var->vp[i];
   }
   var_old->x = var->x - HMIN;
   etc.new_jacobian = 1;
   etc.new_fact = 1; 
   etc.jacobi_status = 1;
   etc.steps_red = 0; 
   etc.h = info->h0;

   set_coeff(&rk); 
   zero_stat();

   etc.gh = info->h0 * rk.a[0][0];

   /* Main loop. */ 
   /*--------------------------------------------------------*/ 
   while (info->xend - var->x > HMIN)
   {
      if ( flag = one_step(&var, &var_old, info, &rk, &etc) )
         return flag;
      if ( info->xend < var->x + etc.h )
      {
         step_fact = (info->xend - var->x)/etc.h; 
         adjust_stepsize(step_fact, &etc);
      }
      ++number_of.steps;
   }

   prob->x = var->x;
   for (i = 0; i < var->n; ++i) 
   {
      prob->v[i] = var->v[i]; 
      prob->vp[i] = var->vp[i];
   }

   work_free(var_old, &etc); 
   return 0;
}
/************ End of sdirk_dae **************************************/


/*---------------------------------------------------------------*/
/* one_step : Takes one step. If appropriate, the stepsize is    */
/*           changed.                                            */
/*---------------------------------------------------------------*/

int one_step(struct problem **var, struct problem **var_old, 
           struct information *info, struct method *rk, 
           struct etcetera *etc)

/*---------------------------------------------------------------*/
/* var     : On entry : Function variables for step n.           */
/*          On exit  : The variables for step n+1                */
/* var_old : On entry : Function variables for step n-1.         */
/*          On exit  : The variables for step n.                 */
/*                                                               */
/* returns : 0 : if one step is taken successfully.              */
/*          1 : The Newton matrix is singular.                   */
/*          2 : The Newton-iterations fails to converge.         */
/*          3 : The local error criteria fails to be satisfied,  */
/*              the step-size becomes too small.                 */
/*          4 : The error tolerance is too stringent.            */
/*---------------------------------------------------------------*/
{
   int i, j, stage, no_success;
   double loc_error, step_fact, roc_fact, sum, error; 
   int n, s1;
   double *h, *gh, *roc, *tol, *v, **vp_stage; 
   struct problem *v_swap;
   extern struct statistic number_of;
/* Set local names on structure variables. */ 
/*-----------------------------------------*/
v = (*var)->v; 
n = (*var)->n; 
s1 = rk->s - 1; 
h = &etc->h; 
gh = &etc->gh; 
roc = &etc->rate_of_convergence;
vp_stage = etc->vp_stage;
tol = etc->tolerance;

if ( compute_tolerance_vector(*var, info, etc) ) 
 return 4;

do
{

 /* If appropriate, the newton matrix is factorized. */ 
 /*                       */
 if (etc ->new_fact) 
 {
  if ( make_newton(*var, etc, info) )
  {
    if ( adjust_stepsize(0.5, etc) )
     return 3;
    if ( make_newton(*var, etc, info) )
     return 1;
  }
 }

 /* Predict starting values for the iterations. */ 
 /*---------------------------------------------*/ 
 prediction(*var, *var_old, rk, etc);

 /* Computes the stage values. */ 
 /*----------------------------*/ 
 no_success = solve_stages(*var, rk, etc);
 roc_fact = ROCK*PES_FAC/(*roc); 
 roc_fact = MAX(0.1, roc_fact); 
 if (no_success)
 {
  ++number_of.notkonv_steps;

  /* The iterations diverge, update the jacobians, or, */ 
  /* if the jacobians is updated, reduce the stepsize. */ 
  /* */ 
  if (etc->jacobi_status)
  {
    etc->new_jacobian = 1;
    etc->new_fact = 1;
  }
  else
  {
    if ( adjust_stepsize( roc_fact, etc) ) 
     return 2;
    etc->steps_red = 2;
  }
 }
   else
 {
     /* The stage values are successfully computed. */ 
     /* Compute the local error.                    */
     /*---------------------------------------------*/ 
     loc_error = 0.0;
     for (i = 0; i < n; ++i)
     {
        sum = 0.0;
        for (stage = 0; stage <= s1; ++stage)
           sum += rk->b_error[stage]*vp_stage[stage][i]; 
        sum *= *h;
        error = fabs(sum)/tol[i];
        loc_error = MAX(loc_error, error);
     }

     /* Assumes 0.5*error from the local error     */ 
     /* and 0.5*error from iterations.             */
     /*--------------------------------------------*/ 
     loc_error *= 2.0;

     if (64.0*loc_error < 1.0) 
        step_fact = PES_FAC*4.0; 
     else
     {
        step_fact = PES_FAC*pow(loc_error, -1.0/3.0);
        step_fact = MAX(0.1, step_fact);
     }
     if (step_fact <= PES_FAC)
     {
        /* The step is rejected. The stepsize is reduced */
        /*                                        */ 
        ++number_of.rejected_steps;
        step_fact = MIN(step_fact, roc_fact);
        if ( adjust_stepsize(step_fact, etc) )
           return 3;
        etc->steps_red = 2;
        no_success = 1;
     }
   }
}
while (no_success);

/* The step is accepted.                                    */
/* Place the old values of the variables in var in var_old. */ 
/*----------------------------------------------------------*/ 
v_swap = *var_old;
*var_old = *var;
*var = v_swap;

/* Calculate the new value for v.                            */
/* Store it in var, together with a new approximation of v', */ 
/* and x(n+1)                                                */
/*-----------------------------------------------------------*/ 
for (i = 0; i < n; ++i)
{
   sum = 0.0;
   for (stage = 0; stage <= s1; ++stage)
     sum += rk->a[s1] [stage]*vp_stage[stage][i]; 
   (*var)->v[i] = v[i] + (*h)*sum;
   (*var)->vp[i] = vp_stage[s1][i];
}
(*var)->x = (*var_old)->x + *h;
  /* If appropriate, change the stepsize.                    */
  /*---------------------------------------------------------*/ 
  step_fact = MIN(step_fact, roc_fact); 
  --etc->steps_red;

  if (etc->h*step_fact > info->hmax)
     step_fact = info->hmax/etc->h;
  if (step_fact > 1.5 && etc->steps_red < 0)
  {
      adjust_stepsize(step_fact, etc);
  }
  else if (step_fact < 1)
  {
     if ( adjust_stepsize(step_fact, etc) )
       return 3;
  } 
  else
  {
     etc->new_fact = 0;
     etc->new_jacobian = 0;
  }
  /* The jacobians are not updated. */ 
  /*--------------------------------*/  
  etc->jacobi_status = 1;

  return 0;
}
/************ End of one_step ***************************************/


/*------------------------------------------------------------------*/ 
/* solve_stages : Solves the Runge-Kutta equations                  */
/*                f(Vi',v + h*a(i,1)*V1,+...+h*a(i,i)*Vi',x+ci*h)=0 */
/*             for the stage values Vi' .....Vi'.                   */
/*------------------------------------------------------------------*/ 

int solve_stages(struct problem *var, struct method *rk, 
 struct etcetera *etc)

/*------------------------------------------------------------------*/ 
/* returns : 0 if the iterations converges successfully.            */
/*           1 if not.                                              */
/*------------------------------------------------------------------*/ 
{
  int i, j, stg;
  double xi, roc;
  double max_roc = 0;
  int s, n;
  double *b, gh, **vp_stage, h;
  double *v, *vp, x;
  extern struct statistic number_of;
  /* Set local names on structure variables. */ 
  /*-----------------------------------------*/
  n = var->n; 
  v = var->v;
  vp = var->vp; 
  x = var->x;

  b = etc->b; 
  gh = etc->gh;
  vp_stage = etc->vp_stage;
  h = etc->h;

  /* Main loop. */ 
  /*------------*/ 
  for (stg = 0; stg < rk->s; ++stg)
  {

     /* Computes b = v + h*a(i,1)*V(1) +... h*a(i,i-1)*V(i-1). */ 
     /*--------------------------------------------------------*/ 
     for (i = 0; i < n; ++i)
     {
        b[i] = v[i];
        for (j = 0; j < stg; ++j)
           b[i] += h * rk->a[stg][j] * vp_stage[j][i];
     }
     xi = x + rk->c[stg] * h;

     /* Solves the equation. Return if not successfully. */
     /*--------------------------------------------------*/
     if ( solve_nonlinear_equations(n, vp_stage[stg], xi, &roc, etc) ) 
     {
        etc->rate_of_convergence = roc;
        return 1;
     }
     max_roc = MAX(max_roc, roc);
  }
   etc->rate_of_convergence = max_roc;
  return 0;
}
/************ End of solve_stages ***********************************/


/*--------------------------------------------------------------*/
/* solve_nonlinear_equation : solves the nonlinear equations    */
/*                             f(vp, gh*vp + b, xi)             */
/*                          with a modified Newton method.      */
/*--------------------------------------------------------------*/

int solve_nonlinear_equations(int n, double *vp, double xi,
                            double *roc, struct etcetera *etc)

/*--------------------------------------------------------------*/
/* n : The dimension of the system.                             */
/* vp : The dependent variable.                                 */
/*      On entry : Some starting values for the iterations.     */
/*      On exit : If successfull, the solution of the eqautions. */
/* xi : xi = xn + c(i)*h.                                       */
/*                                                              */
/* returns : 0 if the iterations converges successfully.        */
/*          1 if the iterations diverges or converges too slow. */
/*--------------------------------------------------------------*/ 
{
   int first_time = 1, near_conv = 0;
   int i, j;
   double error, kappa, last_error, res_error, err_i;
   double *res, *b, *v, *tol, mintol, gh, h;
   extern struct statistic number_of;
 /* Set local names of structure variables. */
 /*-----------------------------------------*/

 gh = etc->gh; 
 h = etc->h;
 b = etc->b;
 mintol = etc->min_tolerance;
 tol =etc->tolerance;
 res = etc->temp1;
 v = etc->temp2;

 /* The iteration loop. */
 /*---------------------*/
 for (*roc = 0.05; (*roc < 1.0) && (*roc < ROCK || near_conv < 2); )
 {
   /* Computes the residual error. */
   /*------------------------------*/ 
   for (i = 0; i < n; ++i)
     v[i] = b[i] + gh*vp[i];
   func(vp, v, xi, res);
   ++number_of.function_evaluations;
  
   /* Solves the linear system. */
   /*---------------------------*/
   backsub(n, etc->newton, etc->pivot, res); 
   ++number_of.solutions;

   /* vp is updated.  */
   /* Computes the displacement error. */ 
   /*----------------------------------*/ 
   error = 0.0;
   for (i = 0; i < n; ++i)
   {
     vp[i] -= res[i];
     err_i = fabs(res[i])*h/tol[i]; 
     error = MAX(error, err_i);
   }

   /* Computes rate of convergence. */
   /*-------------------------------*/ 
   if (first_time)
   {
     first_time = 0;
     if (error < KAPPA)
      return 0;
   } 
   else
   {
     *roc = MAX(0.05,error/last_error);
     if ( fabs(*roc/(1.0 - *roc))*error < KAPPA)
      return 0;
     if (*roc >= ROCK)
      ++near_conv;
   }
   last_error = error;
 }
 return 1;
}
/************ End of solve_nonlinear_equations **********************/

/*------------------------------------------------------------------*/
/* make_newton : Computes the Newton-matrix df/dv, + g*h*df/dv,     */
/*               and performs the LU-decomposition of the matrix.   */
/*------------------------------------------------------------------*/

int make_newton(struct problem *var,
                struct etcetera *etc, struct information *info)
/*------------------------------------------------------------------*/
/* returns: 0 if the factorisation was perfomed successfully.       */
/*          1 if the Newton matrix is singular.                     */
/*------------------------------------------------------------------*/

{
    int i, j;
    int *jacobi_status, *pivot, num_jacobian, new_jacobian, n;
    double *dfdvp, *dfdv, *newton, gh; 
    extern struct statistic number_of;

    /* Set local names of structure variables. */ 
    /*-----------------------------------------*/ 
    n = var->n;
    newton = etc->newton;
    pivot = etc->pivot;
    gh = etc->gh;
    dfdvp = etc->dfdvp;
    dfdv = etc->dfdv;
    jacobi_status = &(etc->jacobi_status); 
    new_jacobian = etc->new_jacobian;
    num_jacobian = info->num_jacobian;

    /* Computes new jacobian if necessary. */ 
    /*-------------------------------------*/ 
    if (new_jacobian)
    {  if (num_jacobian)
          jacobi(var, etc);
       else
          exact_jacobi(var->v, var->vp, var->x, dfdvp, dfdv); 
       *jacobi_status = 0;
       ++number_of.jacobi_evaluations;
    }

    /* Construct the Newton-matrix. */ 
    /*------------------------------*/ 
    for (i = 0; i < n; ++i)
       for (j = 0; j < n; ++j)
          newton[n*i+j] = dfdvp[n*i+j] + gh*dfdv[n*i+j];

    /* Performs the LU-decomposition of the Newton-matrix. */ 
    /*-----------------------------------------------------*/ 
    ++number_of.factorisations;
    if (decomp(n, newton, pivot))
       return 1;
    else
    {  etc->newton_status = 0;
       return 0;
    }
}
/************ End of make_newton ************************************/

/*----------------------------------------------------------*/
/* jacobi : Computes an aproximation to the jacobi-matrices */
/*         df/dv and df/dvp.                                */
/*----------------------------------------------------------*/

int jacobi(struct problem *var, struct etcetera *etc) 
{
   int i, j;
   double delta, absv;
   int n;
   double *v, *vp, x;
   double *res, *res_delta, *dfdv, *dfdvp;
   extern struct statistic nuMber_of;

   /* Set local names of structure variables. */ 
   /*-----------------------------------------*/
   n = var->n; 
   v = var->v; 
   vp = var->vp; 
   x = var->x; 
   dfdv = etc->dfdv;
   dfdvp = etc->dfdvp;
   res = etc->temp1;
   res_delta = etc->temp2;

   /* Start of the computation. */ 
   /*---------------------------*/ 
   func(vp, v, x, res);
   ++number_of.function_evaluations;

   for (j = 0; j < n; ++j)
   {
      /* Computation of df/dvp(j). */
      /*---------------------------*/
      absv = fabs(vp[j]);
      delta = (MAX(1.0, absv))*EPS;
      vp[j] += delta;
      func(vp, v, x, res_delta);
      ++number_of.function_evaluations;
      for (i = 0; i < n; ++i)
         dfdvp[n*i+j] = (res_delta[i] - res[i])/delta;
      vp[j] -= delta;

      /* Computation of df/dv(j). */ 
      /*--------------------------*/ 
      absv = fabs(v[j]);
      delta = (MAX(1.0, absv))*EPS; 
      v[j] += delta;
      func(vp, v, x, res_delta); 
      ++number_of.function_evaluations; 
      for (i = 0; i < n; ++i)
         dfdv[n*i+j] = (res_delta[i] - res[i])/delta;
      v[j] -= delta;
    }
   return 0;
}
/************ End of Jacobi *****************************************/

/*----------------------------------------------------------------*/
/* compute_tolerance_vector : Set tolerance = abserr + relerr*vn. */ 
/*----------------------------------------------------------------*/

int compute_tolerance_vector(struct problem *var,
                           struct information *info, 
                           struct etcetera *etc)

/*------------------------------------------------------------------*/
/* returns : 0 : All elements in the vector are sufficiently large. */ 
/*           1 : The tolerance is near 0 for at least one variables.*/
/*------------------------------------------------------------------*/ 
{
   int i, n;
   double mintol = 1000.0;
   double *tol, *abserr, relerr, *v;

   /* Set local names on structure variables. */ 
   /*-----------------------------------------*/

   n = var->n;
   v = var->v;
   abserr = info->abserr; 
   relerr = info->relerr; 
   tol = etc->tolerance;

   for (i = 0; i < n; ++i)
   {
     tol[i] = abserr[i] + relerr*fabs(v[i]); 
     mintol = MIN(mintol, tol[i]);
   }
   if (mintol < EPS)
     return 1;
   else
     etc->min_tolerance = mintol;

   return 0;
}
/************ End of compute_tolerance_vector ***********************/


/*--------------------------------------------------------------*/
/* prediction : Prediction of the stage_values.                 */
/*             Used as starting values for the iterations       */
/* 1.edition. : Vi = vn.                                        */
/*--------------------------------------------------------------*/

int prediction( struct problem *var, struct problem *var_old, 
               struct method *rk, struct etcetera *etc)

/*--------------------------------------------------------------*/
/* var     : Function variables for step n.                     */
/* var_old : Function variables for step n-1.                   */
/*--------------------------------------------------------------*/ 
{
   int i, stage, n, s;
   double **vp_stage, *vp, h;
   double delta;

   /* Set local names to structure variables. */ 
   /*-----------------------------------------*/ 
   n = var->n;
   vp = var->vp;
   s = rk->s;
   vp_stage = etc->vp_stage;
   h = etc->h;
   /* Main loop. */ 
   /*------------*/ 
   for (i = 0; i < n; ++i)
   {
     delta = (var->vp[i] - var_old->vp[i])/(var->x - var_old->x);
     for (stage = 0; stage < s; ++stage)
     vp_stage[stage][i] = delta*rk->c[stage]*h + vp[i];
   }

   return 0;
}
/************ End of predict ****************************************/



/*---------------------------------------------------------------*/
/* adjust_stepsize : Increase (or decrease) the stepsize with a  */
/*                  factor step_fact.                            */
/*---------------------------------------------------------------*/

int adjust_stepsize(double step_fact, struct etcetera *etc)
{
   etc->h *= step_fact; 
   if (etc->h < HMIN) 
     return 1;
   etc->gh *= step_fact; 
   etc->new_fact = 1; 
   etc->new_jacobian = 0; 
   return 0;
}
/************ End of adjust_stepsize ********************************/



/*---------------------------------------------------------------*/
/* decomp: Gaussian Elimination with Partial Pivoting.           */
/*        References : Golub & van Loan "Matrix Computations",   */
/*        Algorithm 4.4-2.                                       */
/*---------------------------------------------------------------*/

int decomp(int n, double *a, int *piv)

/*---------------------------------------------------------------*/
/* n   : The dimension of the matrix.                            */
/* a   : Pointer to an array of dimension n*n.                   */
/*      On entry : The matrix to be factorised.                  */
/*      On exit  : The LU-decomposition of the matrix.           */
/* piv : The pivot array, of dimension n.                        */
/*---------------------------------------------------------------*/
{
   int i, j, k, p, nml; 
   double maxpivot, temp; 
   nml = n-1;

   /* Start of the factorisation. */ 
   /*-----------------------------*/ 
   for (k = 0; k < nml; ++k)
   {
     maxpivot = 0.0;

      /* Find index of max element in column below the diagonal. */ 
      /*---------------------------------------------------------*/ 
     for (p = k; p < n; ++p)
        if (maxpivot < fabs(a[n*p+k]))
        {
           maxpivot = fabs(a[n*p+k]);
           piv[k] = p;
        }
     /* Test for singularity. */ 
     /*-----------------------*/ 
     if (maxpivot < EPS)
        return 1;

     /* If necessary, swap a(k,j) and a(piv[k],j) for j = 0,..,n-1 */ 
     /*------------------------------------------------------------*/ 
     if (piv[k] != k)
        for (j = 0; j < n; ++j)
        {
           temp = a[n*piv[k]+j]; 
           a[n*piv[k]+j] = a[n*k+j]; 
           a[n*k+j] = temp;
        }

     /* Column elimination. */ 
     /*---------------------*/ 
     for (i = k+1; i < n; ++i) 
     {  a[n*i+k] /= a[n*k+k];
        for (j = k+1; j < n; ++j)
           a[n*i+j] -= a[n*i+k] * a[n*k+j];
     }
  }

  if (fabs(a[n*n-1]) < EPS) return 1;

  /* Everything seems to work, return. */ 
  /*-----------------------------------*/ 
  return 0;
}
/************ End of decomp *****************************************/


/*--------------------------------------------------------------*/
/* backsub: solves the linear equation a*x = b, where a is the  */
/*         factorised matrix, given from the function decomp.   */
/*--------------------------------------------------------------*/

int backsub(int n, double *a, int *piv, double *b)

/*--------------------------------------------------------------*/
/* it  : The dimension of the system.                           */
/* a   : Pointer to an array of dimension n*n.                  */
/*      On entry : The LU-decomposition of the matrix.          */
/* piv : The pivot array of dimension n.                        */
/* b   : Pointer to an array of dimension n.                    */
/*      On entry : The vector b on the right hand side.         */
/*      On exit  : The solution of the equation.                */
/*--------------------------------------------------------------*/
{ int i, j, nm1; 
  double temp; 
  nm1 = n-1;

  /* Swap b[i] and b[piv[i]] */ 
  /*-------------------------*/ 
  for (i = 0; i < nm1; ++i)
     if (piv[i] != i)
     {
        temp = b[piv[i]]; 
        b[piv[i]] = b[i]; 
        b[i] = temp;
     }
   /* Forward elimination. Ly = b. */ 
   /*------------------------------*/ 
   for (i = 0; i < n; ++i)
     for (j = 0; j < i; ++j)
        b[i] -= a[n*i+j]*b[j];

   /* Back substitution. Ux = y*/
   /*--------------------------*/ 
   for (i = nm1; i >= 0; --i)
   {
     for (j = i+1; j < n; ++j) 
        b[i] -= a[n*i+j]*b[j]; 
     b[i] /= a[n*i+i];
   }
   return 0;
}
/************ End of backsub ****************************************/


/*--------------------------------------------------------------*/
/* set_coeff : Set the method coefficients.                     */
/*--------------------------------------------------------------*/

int set_coeff(struct method *rk)
{
   rk ->s = 4;

   rk->a[0][0]= 1.0/4.0;
   rk->a[1][0]= 1.0/7.0;
   rk->a[1][1]= 1.0/4.0;
   rk->a[2][0]= 61.0/144.0;
   rk->a[2][1]= -49.0/144.0;
   rk->a[2][2]= 1.0/4.0;
   rk->a[3][0]= 0.0;
   rk->a[3][1]= 0.0;
   rk->a[3][2]= 3.0/4.0;
   rk->a[3][3]= 1.0/4.0;

   rk->c[0]= 1.0/4.0;
   rk->c[1]= 11.0/28.0;
   rk->c[2]= 1.0/3.0;
   rk->c[3]= 1.0;

   rk->b_error[0]= -61.0/704.0;
   rk->b_error[1]= 49.0/704.0;
   rk->b_error[2]= 3.0/88.0;
   rk->b_error[3]= -3.0/176.0;

   return 0;
}
/************ End of sat_cssff **************************************/

/*--------------------------------------------------------------*/
/* zero_stat : Set all the statistic variables to zero.         */
/*--------------------------------------------------------------*/

int zero_stat()
{
    extern struct statistic number_of; 
    number_of.steps = 0;
    number_of.rejected_steps = 0; 
    number_of.notkonv_steps = 0; 
    number_of.function_evaluations = 0; 
    number_of.jacobi_evaluations = 0; 
    number_of.factorisations = 0; 
    number_of.solutions = 0;

    return 0;
}
/************ End of zero_stat **************************************/


/*--------------------------------------------------------------*/
/* work_alloc : Allocates memory for the work structures.       */
/*--------------------------------------------------------------*/

int work_alloc(int n, struct problem *var, struct etcetera *etc) 
{
    int i;

    /* Allocates memory for the arrays in the structure var. */ 
    /*-------------------------------------------------------*/
    var->v = (double *) malloc( n*sizeof(double) );
    if (!var->v)
    {
      printf("Error in work_alloc: Out of memory: \n");
      exit(1);
    }
    var->vp = (double *) malloc( n*sizeof(double) );
    if (!var->vp)
    {
      printf("Error in work_alloc: Out of memory\n");
      exit (1);
    }

    /* Allocates memory for the arrays in the structure etc. */ 
    /*-------------------------------------------------------*/
    etc->dfdv = (double *) malloc( (n*n)*sizeof(double) );
    if (!etc->dfdv)
    {
      printf("Error in work_alloc: Out of memory\n");
      exit(1);
    }
    etc->dfdvp = (double *) malloc( (n*n)*sizeof(double) );
    if (!etc->dfdvp)
    {
      printf("Error in work_alloc: Out of memory\n");
      exit (1);
    }
    etc->newton = (double *) malloc( (n*n)*sizeof(double) );
    if (!etc->newton)
    {
      printf("Error in work_alloc: Out of memory\n");
      exit (1);
    }
  etc->pivot = (int *) malloc( n*sizeof(int) );
  if (!etc->pivot)
  {
   printf("Error in work_alloc: Out of memory\n");
   exit (1);
  }
  etc->tolerance = (double *) malloc( n*sizeof(double) );
  if (!etc->tolerance)
  {
   printf("Error in work_alloc: Out of memory\n"); 
   exit (1);
  }
  etc->b = (double *) malloc( n*sizeof(double) ); 
  if (!etc->b)
  {
   printf("Error in work_alloc: Out of memory\n"); 
   exit(1);
  }
  etc->temp1 = (double *) malloc( n*sizeof(double) ); 
  if (!etc->temp1)
  {
   printf("Error in work_alloc: Out of memory\n"); 
   exit (1);
  }
  etc->temp2 = (double *) malloc( n*sizeof(double) );
  if (!etc->temp2)
  {
   printf("Error in work_alloc: Out of memory\n");
   exit (1);
  }
  for (i = 0; i < 4; ++i)
  {
   etc->vp_stage[i] = (double *) malloc( n*sizeof(double) );
   if (!etc->vp_stage[i])
   {
    printf("Error in work_alloc: Out of memory\n");
    exit (1);
   }
  }
  return 0;
}
/************ End of nork_allno *************************************/


/*--------------------------------------------------------------*/
/* work_free : Free the memory used for the work structures.    */
/*--------------------------------------------------------------*/

int work_free(struct problem *var, struct etcetera *etc)
{
  int i;
/* free(var->v);
  free(var->vp); */

  free(etc->dfdv); 
  free(etc->dfdvp); 
  free(etc->newton); 
  free(etc->pivot); 
  free(etc->tolerance);
  free(etc->b);
  free(etc->temp1); 
  free(etc->temp2); 
  for (i = 0; i < 4; ++i)
   free(etc->vp_stage[i]);
  return 0;
}
/************ End of work_free **********************************/
