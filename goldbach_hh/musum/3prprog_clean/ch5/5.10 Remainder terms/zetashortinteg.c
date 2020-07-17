/*
Compile with

gcc -o zetashortinteg zetashortinteg.c -larb -lflint
*/

/*
  Rigorous numerical integration (with fast convergence for piecewise
  holomorphic functions) using Gauss-Legendre quadrature and adaptive
  subdivision.

  Author: Fredrik Johansson.
  Modified by H. A. Helfgott
  This file is in the public domain.

  Todo:
  * Allow relative instead of absolute tolerance (based on the computed
    values of the integrand)
  * Multithreaded evaluation
  * More examples / test integrals
  * Improvements to the adaptive evaluation and error bounding strategy
  * Better work limits / abort criteria
  * Other quadrature formulas (Clenshaw-Curtis, double exponential...)
  * When this is more stable, make the relevant parts library routines
*/

#include <string.h>
#include "flint/profiler.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"


/* ------------------------------------------------------------------------- */
/*  Integrands                                                               */
/* ------------------------------------------------------------------------- */

/*
  A function to evaluate the integrand res = f(z). It can be assumed that
  res and z are not aliased.

  An integrand will get called in two ways. If holomorphic = 0, the function
  is evaluated normally on the integration path (either at a single point
  or on a subinterval). In this case, f is treated as a pointwise defined
  function and can have arbitrary discontinuities.

  If holomorphic = 1, the function is evaluated on a domain surrounding
  a subinterval of the integration path (for the purpose of computing an error
  bound for a quadrature formula). In this case, the evaluation function MUST
  verify that f is holomorphic on this domain and set res to a non-finite
  value otherwise.

  If f is built from field operations and meromorphic functions, then
  no special action is necessary since division by zero automatically will
  produce an infinite enclosure.

  However, manual action is needed for bounded functions with branch cuts.
  For example, when evaluating res = sqrt(z) with holomorphic = 1, then
  res must be set to an non-finite value if z overlaps with the standard
  branch cut [-inf,0].

  The easiest solution is to use holomorphy-testing versions of atomic
  functions. As examples, we define sqrt() and piecewise holomorphic
  extensions of the piecewise real analytic functions abs() and floor()
  below. See the example integrals for usage.
*/

typedef void (*integrand_t)
    (acb_t res, void * param, const acb_t z, int holomorphic, slong prec);


/* ------------------------------------------------------------------------- */
/*  Integrands                                                       */
/* ------------------------------------------------------------------------- */


void fac(acb_t res, slong p, const acb_t s, slong prec)
/* sets res to 1 - p^{-s} */
{
  acb_t pc, ms;

  acb_init(pc);   acb_init(ms);

  acb_set_si(pc,p); acb_neg(ms,s); 
  acb_pow(res,pc,ms,prec);
	  
  acb_sub_si(res,res,1,prec);
  acb_neg(res,res);

  acb_clear(pc);   acb_clear(ms);
}

void fillisprime(short int *isprime, long N)
/* sets isprime[n]=0 if n is composite, isprime[n]=1 i n is prime
   for 2<=n<N */
{
  long i,j;

  for(i=2; i<N; i++)
    isprime[i]=1;
  for(i=2; i*i<N; i++)
    if(isprime[i])
      for(j=i; j*i<N; j++)
        isprime[j*i]=0;
}

void P(acb_t res, const acb_t s, short int *isprime, int m, slong prec)
/* sets res to P_m(s) = \prod_{p\leq m} (1-p^{-s}) */
 {
   acb_t facn;
   int p;
   
   if(m<2) {
     acb_set_si(res,(slong) 1);
     return;
   }
   
   acb_init(facn);

   fac(res,2,s,prec);

   for(p=3; p<=m; p+=2)
     if(isprime[p]) {
       fac(facn,p,s,prec); acb_mul(res,res,facn,prec);
     }

   acb_clear(facn);
 }

short int *isprime;

void
f_myintegr(acb_t res, void * param, const acb_t z, int holomorphic, slong prec)
/* sets res to |P_19(z+1/2) P_19(z+1) zeta(z+1/2) zeta(z+1)|/|z|^2 
 */
{
  acb_t sr,s1,s2, half;

  acb_init(sr); acb_init(s1); acb_init(s2); acb_init(half);

  acb_set_d(half,0.5);
  acb_add(s1,z,half,prec);
  acb_add_si(s2,z,1,prec);
  
  acb_zeta(res, s1, prec);

  P(sr,s1,isprime,19,prec);
  acb_mul(res,res,sr,prec);
  
  acb_zeta(sr, s2, prec);
  acb_mul(res,res,sr,prec);

  P(sr,s2,isprime,19,prec);
  acb_mul(res,res,sr,prec);
  
  acb_mul(sr,z,z,prec);
  acb_div(res,res,sr,prec);

  acb_abs(acb_realref(res), res, prec);
  arb_zero(acb_imagref(res));
  
  acb_clear(sr); acb_clear(s1); acb_clear(s2); acb_clear(half);
}

void
f_myintegr2(acb_t res, void * param, const acb_t z, int holomorphic, slong prec)
/* sets res to |P_{11}(z+1/2) P_{11}(z+1/4) zeta(z+1/2) zeta(z+1/4)|/|z|^2 
 */  
{
  acb_t sr,s1,s2, half, quart;

  acb_init(sr); acb_init(s1); acb_init(s2); acb_init(half); acb_init(quart);

  acb_set_d(half,0.5);   acb_set_d(quart,0.25);
  acb_add(s1,z,half,prec);
  acb_add(s2,z,quart,prec);
  
  acb_zeta(res, s1, prec);

  P(sr,s1,isprime,11,prec);
  acb_mul(res,res,sr,prec);
  
  acb_zeta(sr, s2, prec);
  acb_mul(res,res,sr,prec);

  P(sr,s2,isprime,11,prec);
  acb_mul(res,res,sr,prec);
  
  acb_mul(sr,z,z,prec);
  acb_div(res,res,sr,prec);

  acb_abs(acb_realref(res), res, prec);
  arb_zero(acb_imagref(res));
  
  acb_clear(sr); acb_clear(s1); acb_clear(s2); acb_clear(half); acb_clear(quart);
}



/* ------------------------------------------------------------------------- */
/*  Quadrature nodes                                                         */
/* ------------------------------------------------------------------------- */

/*
  Gauss-Legendre quadrature nodes are cached to speed up multiple integrations
  and adaptive subdivision. The steps of 2^(n/2) used here give slightly better
  performance than steps of 2^n (we use at most 1.4x more points than
  needed and not 2x more points) but may require more precomputation.
*/

#define GL_STEPS 38

const slong gl_steps[GL_STEPS] = {1, 2, 4, 6, 8, 12, 16, 22, 32, 46, 64,
    90, 128, 182, 256, 362, 512, 724, 1024, 1448, 2048, 2896, 4096,
    5792, 8192, 11586, 16384, 23170, 32768, 46340, 65536, 92682,
    131072, 185364, 262144, 370728, 524288, 741456};

slong gl_prec[GL_STEPS] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

arb_ptr gl_nodes[GL_STEPS];
arb_ptr gl_weights[GL_STEPS];

void gl_cleanup()
{
    slong i;
    for (i = 0; i < GL_STEPS; i++)
    {
        if (gl_prec[i] != 0)
        {
            _arb_vec_clear(gl_nodes[i], (gl_steps[i] + 1) / 2);
            _arb_vec_clear(gl_weights[i], (gl_steps[i] + 1) / 2);
        }
    }
}

/* Compute GL node and weight of index k for n = gl_steps[i]. Cached. */
void gl_node(arb_t x, arb_t w, slong i, slong k, slong prec)
{
    slong n, kk, jj, wp;

    if (i < 0 || i >= GL_STEPS || prec < 2)
        flint_abort();

    n = gl_steps[i];

    if (k < 0 || k >= n)
        flint_abort();

    if (2 * k < n)
        kk = k;
    else
        kk = n - 1 - k;

    if (gl_prec[i] < prec)
    {
        if (gl_prec[i] == 0)
        {
            gl_nodes[i] = _arb_vec_init((n + 1) / 2);
            gl_weights[i] = _arb_vec_init((n + 1) / 2);
        }

        wp = FLINT_MAX(prec, gl_prec[i] * 2 + 30);

        for (jj = 0; 2 * jj < n; jj++)
        {
            arb_hypgeom_legendre_p_ui_root(gl_nodes[i] + jj,
                gl_weights[i] + jj, n, jj, wp);
        }

        gl_prec[i] = wp;
    }

    if (2 * k < n)
        arb_set_round(x, gl_nodes[i] + kk, prec);
    else
        arb_neg_round(x, gl_nodes[i] + kk, prec);

    arb_set_round(w, gl_weights[i] + kk, prec);
}

/* ------------------------------------------------------------------------- */
/*  Adaptive quadrature                                                      */
/* ------------------------------------------------------------------------- */

/*
Attempts to compute integral on [a,b] without subdivisions,
using at most O(prec) evaluations.

A return value of 0 indicates that [a,b] should be subdivided.

A return value of 1 indicates that the integral was computed successfully:
not necessarily to the specified tolerance (tol), but further subdivisions
are unlikely to give a much better result.

Strategy: we first evaluate f([a,b]) and check if the trivial enclosure
is good enough. If not, we try to find a good Gauss-Legendre rule.
For the interval [-1,1], the error of the n-point GL rule is bounded by

    (64/15) M / ((rho-1) rho^(2n-1))

where f(z) is analytic inside the ellipse with foci +/- 1 and semiaxes X, Y
such that 1 < rho = X + Y = X + sqrt(X^2 - 1) and |f(z)| <= M on this ellipse
(see Trefethen, "Is Gauss Quadrature Better than Clenshaw-Curtis?").

For an arbitrary interval, we use int_a^b f(x) = int_-1^1 g(x)
where g(x) = (b-a)/2  f((b-a)/2 x + (a+b)/2). With I = [+/-X] + [+/-Y]i,
this means we evaluate (b-a)/2  f((b-a)/2 I + (a+b)/2) to get the bound M
(todo: reduce the wrapping effect of rotating the ellipse when the path
is not rectilinear).

We search for an X that makes the error small by trying steps 2^(2^k).
Larger X will give smaller 1/rho^(2n-1) but larger M. If we try
successive larger values of k, we can abort when M = inf since this either
means that we have hit a singularity or a branch cut or that overestimation
in the evaluation of f is becoming too severe.
*/

int
quad_auto_deg(acb_t res, integrand_t f, void * param,
        const acb_t a, const acb_t b, const mag_t tol, slong prec)
{
    acb_t mid, delta, wide;
    mag_t tmpm;
    int success, real_err;

    success = 0;

    acb_init(mid);
    acb_init(delta);
    acb_init(wide);
    mag_init(tmpm);

    /* delta = (b-a)/2 */
    acb_sub(delta, b, a, prec);
    acb_mul_2exp_si(delta, delta, -1);

    /* mid = (a+b)/2 */
    acb_add(mid, a, b, prec);
    acb_mul_2exp_si(mid, mid, -1);

    /* wide = mid +- [delta] */
    acb_set(wide, mid);
    arb_get_mag(tmpm, acb_realref(delta));
    arb_add_error_mag(acb_realref(wide), tmpm);
    arb_get_mag(tmpm, acb_imagref(delta));
    arb_add_error_mag(acb_imagref(wide), tmpm);

    /* Direct evaluation: integral = (b-a) * f([a,b]). */
    f(res, param, wide, 0, prec);
    acb_mul(res, res, delta, prec);
    acb_mul_2exp_si(res, res, 1);

    mag_hypot(tmpm, arb_radref(acb_realref(res)), arb_radref(acb_imagref(res)));
    real_err = (acb_is_finite(res) && acb_is_real(res));

    if (acb_contains_zero(delta) || mag_cmp(tmpm, tol) < 0)
    {
        success = 1;
    }
    else
    {
        acb_t s, v;
        mag_t M, X, Y, RHO, err, t;
        slong k, Xexp;
        slong i, n, best_n;

        acb_init(s);
        acb_init(v);
        mag_init(M);
        mag_init(X);
        mag_init(Y);
        mag_init(RHO);
        mag_init(t);
        mag_init(err);

        best_n = -1;

        mag_inf(err);

        for (Xexp = 0; Xexp < prec; Xexp += FLINT_MAX(1, Xexp))
        {
            mag_one(X);
            mag_mul_2exp_si(X, X, Xexp + 1);

            /* RHO = X + sqrt(X^2 - 1)  (lower bound) */
            mag_mul_lower(RHO, X, X);
            mag_one(t);
            mag_sub_lower(RHO, RHO, t);
            mag_sqrt_lower(RHO, RHO);
            mag_add_lower(RHO, RHO, X);

            /* Y = sqrt(X^2 - 1)  (upper bound) */
            mag_mul(Y, X, X);
            mag_one(t);
            mag_sub(Y, Y, t);
            mag_sqrt(Y, Y);

            acb_zero(wide);
            mag_set(arb_radref(acb_realref(wide)), X);
            mag_set(arb_radref(acb_imagref(wide)), Y);

            /* transform to [a,b] */
            acb_mul(wide, wide, delta, prec);
            acb_add(wide, wide, mid, prec);

            f(v, param, wide, 1, prec);

            /* no chance */
            if (!acb_is_finite(v))
                break;

            /* M = (b-a)/2  |f| */
            acb_get_mag(M, v);
            acb_get_mag(tmpm, delta);
            mag_mul(M, M, tmpm);

            /* Search for the smallest n that gives err < tol (if possible) */
            for (i = 0; i < GL_STEPS && gl_steps[i] < 0.5 * prec + 10; i++)
            {
                n = gl_steps[i];

                /* (64/15) M / ((rho-1) rho^(2n-1)) */
                mag_pow_ui_lower(t, RHO, 2 * n - 1);
                mag_one(tmpm);
                mag_sub_lower(tmpm, RHO, tmpm);
                mag_mul_lower(t, t, tmpm);
                mag_mul_ui_lower(t, t, 15);
                mag_div(t, M, t);
                mag_mul_2exp_si(t, t, 6);

                if (mag_cmp(t, tol) < 0)
                {
                    success = 1;

                    /* The best so far. */
                    if (best_n == -1 || n < best_n)
                    {
                        mag_set(err, t);
                        best_n = n;
                    }

                    /* Best possible n. */
                    if (n == 1)
                        break;
                }
            }
        }

        /* Evaluate best found Gauss-Legendre quadrature rule. */
        if (success)
        {
            arb_t x, w;
            arb_init(x);
            arb_init(w);

#if 0
            flint_printf("  {deg %ld for ", best_n);
            acb_printn(a, 20, ARB_STR_NO_RADIUS); flint_printf("   ");
            acb_printn(b, 20, ARB_STR_NO_RADIUS); flint_printf("}\n");
#endif

            if (best_n == -1)
                flint_abort();

            for (i = 0; i < GL_STEPS; i++)
                if (gl_steps[i] == best_n)
                    break;

            acb_zero(s);

            for (k = 0; k < best_n; k++)
            {
                gl_node(x, w, i, k, prec);
                acb_mul_arb(wide, delta, x, prec);
                acb_add(wide, wide, mid, prec);
                f(v, param, wide, 0, prec);
                acb_addmul_arb(s, v, w, prec);
            }

            acb_mul(res, s, delta, prec);

            if (real_err)
                arb_add_error_mag(acb_realref(res), err);
            else
                acb_add_error_mag(res, err);

            arb_clear(x);
            arb_clear(w);
        }

        acb_clear(s);
        acb_clear(v);
        mag_clear(M);
        mag_clear(X);
        mag_clear(Y);
        mag_clear(RHO);
        mag_clear(t);
        mag_clear(err);
    }

    acb_clear(mid);
    acb_clear(delta);
    acb_clear(wide);
    mag_clear(tmpm);

    return success;
}

/* Adaptive subdivision on [a,b]. */
void
quad_subdiv(acb_t res, integrand_t f, void * param,
    const acb_t a, const acb_t b, const mag_t tol, slong prec)
{
    acb_ptr as, bs;
    acb_t s, t;
    slong maxdepth, num, eval, found, num_most;

    acb_init(s);
    acb_init(t);

    maxdepth = 10 * prec;
    maxdepth = FLINT_MAX(maxdepth, 1);

    /* todo: allocate dynamically */
    as = _acb_vec_init(maxdepth);
    bs = _acb_vec_init(maxdepth);
    acb_set(as, a);
    acb_set(bs, b);

    num = 1;
    num_most = 1;
    eval = found = 0;

    acb_zero(s);

    while (num >= 1)
    {
        if (quad_auto_deg(t, f, param, as + num - 1, bs + num - 1, tol, prec) == 1)
        {
            acb_add(s, s, t, prec);
            eval++;
#if 0
            flint_printf("ok [%ld]: ", eval);
            acb_printn(as + num - 1, 10, ARB_STR_NO_RADIUS);
            flint_printf("   ");
            acb_printn(bs + num - 1, 10, ARB_STR_NO_RADIUS);
            flint_printf("\n");
#endif
            num--;
        }
        else
        {
            if (num == maxdepth - 1)
            {
                flint_printf("too many subdivisions!\n");
                flint_abort();
            }

            /* add [a,mid(x)], [mid(x),b] */
            /* interval [num] becomes [mid, b] */
            acb_set(bs + num, bs + num - 1);
            acb_add(as + num, as + num - 1, bs + num - 1, prec);
            acb_mul_2exp_si(as + num, as + num, -1);
            /* interval [num-1] becomes [a, mid] */
            acb_set(bs + num - 1, as + num);
            num++;
            num_most = FLINT_MAX(num, num_most);
        }
    }

    acb_set(res, s);

    _acb_vec_clear(as, maxdepth);
    _acb_vec_clear(bs, maxdepth);
    acb_clear(s);
    acb_clear(t);
}

/* ------------------------------------------------------------------------- */
/*  Main test program                                                        */
/* ------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
    acb_t s, t, a, b;
    mag_t tol;
    slong prec;
    int integral;
    int i, twice;

    flint_printf("Compute integrals using subdivision and Gauss-Legendre quadrature.\n");
    flint_printf("Usage: quadrature [-prec p] [-tol eps] [-twice]\n\n");
    flint_printf("-prec p    - precision in bits (default p = 333)\n");
    flint_printf("-tol eps   - approximate absolute error goal (default 2^-p)\n");
    flint_printf("-twice     - run twice (to see overhead of computing nodes)\n");
    flint_printf("\n\n");

    prec = 333;
    twice = 0;

    acb_init(a);
    acb_init(b);
    acb_init(s);
    acb_init(t);
    mag_init(tol);

    isprime = (short int *)  calloc(101,sizeof(short int));
    fillisprime(isprime,101);
    
    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-prec"))
        {
            prec = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-twice"))
        {
            twice = 1;
        }
        else if (!strcmp(argv[i], "-tol"))
        {
            arb_t x;
            arb_init(x);
            arb_set_str(x, argv[i+1], 10);
            arb_get_mag(tol, x);
            arb_clear(x);
        }
    }

    if (mag_is_zero(tol))
        mag_set_ui_2exp_si(tol, 1, -prec);

    for (integral = 0; integral <= 1; integral++)
    {
        for (i = 0; i < 1 + twice; i++)
        {
            TIMEIT_ONCE_START
            switch (integral)
            {
	    case 0:
                flint_printf("Computing I = 1/(2 pi i) int |P_{19}(s+1/2) P_{19}(s+1)| |zeta(s+1/2) zeta(s+1)/s^2| ds\n\t(straight path from -8i - 1/4 to -8i - 1/4)\n");
                acb_zero(s);

		acb_set_d_d(a, -0.25, -8.0);
                acb_set_d_d(b, -0.25, 0.0);
                quad_subdiv(t, f_myintegr, NULL, a, b, tol, prec);
		acb_abs(acb_realref(t), t, prec);
		arb_zero(acb_imagref(t));
                acb_add(s, s, t, prec);
		
                acb_const_pi(t, prec);
                acb_div(s, s, t, prec);
                break;

	    case 1:
                flint_printf("Computing I = 1/(2 pi i) int |P_{11}(s+1/2) P_{11}(s+1)| |zeta(s+1/2) zeta(s+1/4)/s^2| ds\n\t(straight paths from -200i to -0.005 to 200i)\n");
                acb_zero(s);

		acb_set_d_d(a, 0, -200.0);
                acb_set_d_d(b, -0.005, 0.0);
                quad_subdiv(t, f_myintegr2, NULL, a, b, tol, prec);
		acb_abs(acb_realref(t), t, prec);
		arb_zero(acb_imagref(t));
                acb_add(s, s, t, prec);
		
                acb_const_pi(t, prec);
                acb_div(s, s, t, prec);
                break;

            default:
                abort();
            }
            TIMEIT_ONCE_STOP
        }
        flint_printf("I = ");
        acb_printn(s, 3.333 * prec, 0);
        flint_printf("\n\n");
    }

    acb_clear(a);
    acb_clear(b);
    acb_clear(s);
    acb_clear(t);
    mag_clear(tol);

    gl_cleanup();
    flint_cleanup();
    free(isprime);
    return 0;
}

