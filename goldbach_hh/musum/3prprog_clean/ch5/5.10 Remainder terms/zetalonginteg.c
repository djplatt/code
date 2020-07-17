#include "acb_dirichlet.h"
#include "flint/profiler.h"
#include "acb_calc.h"

/*
Compile with

gcc -o zetalonginteg zetalonginteg.c -larb -lflint
 */

/*
  Rigorous numerical integration on long vertical intervals.

  Author: Fredrik Johansson
  Modified by Harald Helfgott
  This file is in the public domain

 */

void fac(acb_t res, slong p, const acb_t s, slong prec)
/* Bound |F(s)| (crude estimate to bound Taylor tail). */
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

void P(acb_t res, const acb_t s, short int *isprime, int M, slong prec)
/* sets res to P_M(s) = \prod_{p\leq M} (1-p^{-s}) */
 {
   acb_t facn;
   int p;
   
   if(M<2) {
     acb_set_si(res,(slong) 1);
     return;
   }
   
   acb_init(facn);

   fac(res,2,s,prec);

   for(p=3; p<=M; p+=2)
     if(isprime[p]) {
       fac(facn,p,s,prec); acb_mul(res,res,facn,prec);
     }

   acb_clear(facn);
 }

void
bound(mag_t res, const acb_t s, double c1, double c2, short *isprime, int M)
/* gives upper bound on
|zeta(s+c1)| |zeta(s+c2)| |P_M(s+c1)| |P_M(s+c2)|/|s|^2
 */
{
  acb_t t,z;
  mag_t u;

  acb_init(t);     acb_init(z);
  mag_init(u);
  
  acb_set_d(t, c1);
  acb_add(t, t, s, MAG_BITS);
  acb_dirichlet_zeta_bound(res, t);
  
  P(z,t,isprime,M,MAG_BITS);
  acb_get_mag(u,z);
  mag_mul(res,res,u);

  acb_set_d(t, c2);
  acb_add(t, t, s, MAG_BITS);
  acb_dirichlet_zeta_bound(u, t);
  mag_mul(res,res,u);
  
  P(z,t,isprime,M,MAG_BITS);
  acb_get_mag(u,z);
  mag_mul(res,res,u);
  
  acb_get_mag_lower(u, s);
  mag_mul_lower(u, u, u);
  mag_div(res, res, u);
       
    acb_clear(t);     acb_clear(z);
    mag_clear(u);
}

void fac_ser(acb_poly_t res, slong p, const acb_poly_t sser, slong N, slong prec)
{
  acb_t mlp;
  acb_poly_t expo;
  
  acb_init(mlp);  acb_poly_init(expo); 

  acb_set_si(mlp,p);
  acb_log(mlp,mlp,prec);
  acb_neg(mlp,mlp);
  
  acb_poly_scalar_mul(expo,sser,mlp,prec); 
  acb_poly_exp_series(res,expo,N,prec);
	  
  acb_poly_add_si(res,res,-1,prec);
  acb_poly_neg(res,res);

  acb_clear(mlp); acb_poly_clear(expo); 
}

void P_ser(acb_poly_t res, acb_poly_t sser, slong N, short *isprime,
	     int M, slong prec)
{
   acb_poly_t facn;
   int p;
   
   if(M<2) {
     acb_poly_one(res);
     return;
   }
   
   acb_poly_init(facn);

   fac_ser(res,2,sser,N,prec);

   for(p=3; p<=M; p+=2)
     if(isprime[p]) {
       fac_ser(facn,p,sser,N,prec);  acb_poly_mullow(res,res,facn,N,prec);
     }
   
   acb_poly_clear(facn);
}

/* Compute order-N Taylor series of F at the point s. */
void series(acb_poly_t res, const acb_t s, slong N,
	    double c1, double c2, short *isprime, int M,
	    slong prec)
{
  acb_poly_t f, g, h, gres, hres, numo;
  acb_t one, dum;

    acb_poly_init(f);
    acb_poly_init(g);
    acb_poly_init(h);
    acb_poly_init(gres);
    acb_poly_init(hres);     acb_poly_init(numo);
    acb_init(one);
    acb_init(dum);
    
    /* f = s+x */
    acb_poly_set_coeff_acb(f, 0, s);
    acb_poly_set_coeff_si(f, 1, 1);

    /* g = f + c1 */
    acb_set_d(dum, c1);    
    acb_poly_one(g);
    acb_poly_set_coeff_acb(g, 0, dum);
    acb_poly_add(g, g, f, prec);
    
    /* h = f + c2 */
    acb_set_d(dum, c2);    
    acb_poly_one(h);
    acb_poly_set_coeff_acb(h, 0, dum);
    acb_poly_add(h, h, f, prec);

    acb_one(one);
    acb_poly_zeta_series(gres, g, one, 0, N, prec);
    acb_poly_zeta_series(hres, h, one, 0, N, prec);
    acb_poly_mullow(numo, gres, hres, N, prec);

    P_ser(gres, g, N, isprime, M, prec);
    P_ser(hres, h, N, isprime, M, prec);    
    acb_poly_mullow(numo, numo, gres, N, prec);
    acb_poly_mullow(numo, numo, hres, N, prec);
	    
    acb_poly_mullow(f, f, f, N, prec);
    acb_poly_div_series(res, numo, f, N, prec);
    
    acb_poly_clear(f);
    acb_poly_clear(g);  acb_poly_clear(h);
    acb_poly_clear(gres);  acb_poly_clear(hres);
    acb_poly_clear(numo);
    acb_clear(one);
}

/*
Compute Taylor model (P, err) for F(s) such that |P(s-m) - F(s)| <= err for
s = sigma + it with a <= t <= b, m = sigma + (a+b)/2 i.
*/
void taylor(acb_poly_t P, mag_t err,
    const arb_t sigma, const arb_t a, const arb_t b,
	    double c1, double c2, short *isprime, int M, mag_t tol, slong prec)
{
    acb_t m, wide;
    arb_t delta;
    mag_t C, D, r, R, err2, tmpm;
    slong k, N;

    acb_init(m);
    acb_init(wide);
    arb_init(delta);
    mag_init(r);
    mag_init(R);
    mag_init(C);
    mag_init(D);
    mag_init(err2);
    mag_init(tmpm);

    /* r = inner radius, R = outer radius for bound */
    arb_sub(delta, b, a, prec);
    arb_get_mag(r, delta);
    mag_mul_2exp_si(r, r, -1);
    mag_mul_2exp_si(R, r, 1);

    /* m = midpoint of interval */
    arb_set(acb_realref(m), sigma);
    arb_add(acb_imagref(m), a, b, prec);
    arb_mul_2exp_si(acb_imagref(m), acb_imagref(m), -1);

    acb_set(wide, m);
    acb_add_error_mag(wide, R);
    bound(C, wide, c1, c2, isprime, M);

    /* Taylor series error: C D^N / (1 - D),  R = r / R */
    mag_div(D, r, R);

    for (N = 2; ; N++)
    {
        mag_geom_series(err, D, N);
        mag_mul(err, err, C);
        mag_mul_2exp_si(tmpm, tol, -1);
        if (mag_cmp(err, tmpm) < 0)
            break;
    }

    /* Compute Taylor polynomial */
    series(P, m, N, c1, c2, isprime, M, prec);

    /* Economize Taylor polynomial */
    for (k = acb_poly_length(P) - 1; k > 0; k--)
    {
        /* Removing term c*x^k  adds an error |c| r^k  */
        acb_get_mag(err2, P->coeffs + k);
        mag_pow_ui(tmpm, r, k);
        mag_mul(err2, err2, tmpm);

        mag_add(tmpm, err, err2);
        if (mag_cmp(tmpm, tol) < 0)
        {
            acb_poly_set_coeff_si(P, k, 0);
            mag_add(err, err, err2);
        }
    }

    acb_clear(m);
    acb_clear(wide);
    arb_clear(delta);
    mag_clear(r);
    mag_clear(R);
    mag_clear(C);
    mag_clear(D);
    mag_clear(err2);
    mag_clear(tmpm);
}

/* Naive integration of |P(it)|, -(b-a)/2 <= t <= (b-a)/2 */
void
integrate(arb_t I, const acb_poly_t P, mag_t err, const arb_t a, const arb_t b, slong N, slong prec)
{
    acb_t s;
    arb_t t, u, h;
    slong k;

    acb_init(s);
    arb_init(t);
    arb_init(u);
    arb_init(h);

    /* step size h = (b-a)/N */
    arb_sub(h, b, a, prec);
    arb_div_ui(h, h, N, prec);
    arb_mul_2exp_si(h, h, -1);

    for (k = -N; k < N; k++)
    {
        arb_zero(acb_realref(s));
        arb_mul_si(t, h, k, prec);
        arb_mul_si(u, h, k + 1, prec);
        arb_union(acb_imagref(s), t, u, prec);
        acb_poly_evaluate(s, P, s, 53);

        acb_abs(t, s, prec);
        arb_add_error_mag(t, err);

        arb_add(I, I, t, prec);
    }

    arb_abs(t, h);
    arb_mul(I, I, t, prec);

    acb_clear(s);
    arb_clear(t);
    arb_clear(u);
    arb_clear(h);
}

/* Square root function on C with detection of the branch cut. */
void
acb_holomorphic_sqrt(acb_ptr res, const acb_t z, int holomorphic, slong prec)
{
    if (!acb_is_finite(z) || (holomorphic &&
        arb_contains_zero(acb_imagref(z)) &&
        arb_contains_nonpositive(acb_realref(z))))
    {
        acb_indeterminate(res);
    }
    else
    {
        acb_sqrt(res, z, prec);
    }
}

int
f_eval(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t t, u;
    acb_poly_struct * AB;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    AB = (acb_poly_struct *) param;

    acb_init(t);
    acb_init(u);

    acb_poly_evaluate(t, AB + 0, z, prec);
    acb_mul(t, t, t, prec);
    acb_poly_evaluate(u, AB + 1, z, prec);
    acb_mul(u, u, u, prec);
    acb_add(t, t, u, prec);
    acb_holomorphic_sqrt(res, t, order != 0, prec);

    acb_clear(t);
    acb_clear(u);

    return 0;
}

/* Less naive integration of |P(it)|, -(b-a)/2 <= t <= (b-a)/2 */
void
integrate2(arb_t I, const acb_poly_t P, mag_t err, const arb_t a, const arb_t b, const mag_t tol, slong prec)
{
    acb_poly_struct AB[2];
    acb_t c, d, s;
    slong k;
    mag_t err2;

    acb_poly_init(AB + 0);
    acb_poly_init(AB + 1);
    acb_init(c);
    acb_init(d);
    acb_init(s);
    mag_init(err2);

    for (k = 0; k < P->length; k++)
    {
        acb_poly_get_coeff_acb(c, P, k);

        if (k % 4 == 1)
            acb_mul_onei(c, c);
        else if (k % 4 == 2)
            acb_neg(c, c);
        else if (k % 4 == 3)
            acb_div_onei(c, c);

        arb_zero(acb_imagref(d));
        arb_set(acb_realref(d), acb_realref(c));
        acb_poly_set_coeff_acb(AB + 0, k, d);
        arb_set(acb_realref(d), acb_imagref(c));
        acb_poly_set_coeff_acb(AB + 1, k, d);
    }

    acb_zero(c);
    acb_zero(d);
    arb_sub(acb_realref(c), b, a, prec);
    acb_mul_2exp_si(c, c, -1);
    acb_set(d, c);
    acb_neg(c, c);

/*
    acb_poly_printd(AB + 0, 30); printf("\n\n");
    acb_poly_printd(AB + 1, 30); printf("\n\n");
*/

    acb_calc_integrate(s, f_eval, &AB, c, d, prec, tol, NULL, prec);

    /* add approximation error bound: |b-a|*err */
    acb_get_mag(err2, d);
    mag_mul_2exp_si(err2, err2, 1);
    mag_mul(err2, err2, err);
    acb_add_error_mag(s, err2);

    arb_set(I, acb_realref(s));

    acb_poly_clear(AB + 0);
    acb_poly_clear(AB + 1);
    acb_clear(c);
    acb_clear(d);
    acb_clear(s);
    mag_clear(err2);
}

int main(int argc, char *argv[])
{
    acb_poly_t P;
    arb_t sigma, a, b, I, I2, I3, pi;
    mag_t tol, err;
    slong prec, M;
    double k, h;
    double c1, c2;
    short *isprime;
    int m, k0, k1;
    
    acb_poly_init(P);
    arb_init(sigma); arb_init(pi);
    arb_init(a);
    arb_init(b);
    arb_init(I);
    arb_init(I2);
    arb_init(I3);
    mag_init(tol);
    mag_init(err);

    prec = 53;
    M = 32768;
    mag_set_d(tol, 1e-12);
    h = 0.5;

    arb_const_pi(pi, 256);

    if(argc<7) {
      m=19;
      c1 = 1.0; c2 = 0.5;
      arb_set_d(sigma, -0.25);
      k0=8; k1 = 40000;
    } else {
      m = atoi(argv[1]);
      c1 = atof(argv[2]); c2 = atof(argv[3]);
      arb_set_d(sigma, atof(argv[4]));
      k0 = atoi(argv[5]); k1 = atoi(argv[6]);
    }
    
    isprime = (short *) calloc(m+1,sizeof(short));
    fillisprime(isprime,m+1);
		
    for (k = k0; k < k1; k += h)
    {
        arb_set_d(a, k);
        arb_set_d(b, k + h);

        taylor(P, err, sigma, a, b, c1, c2, isprime, m, tol, prec);

        /* integrate(I2, P, err, a, b, M, prec); */
        integrate2(I3, P, err, a, b, tol, prec);

        arb_add(I, I, I3, prec);
        printf("[%f, %f] ", k, k + h);
        /*
        arb_printn(I2, 10, 0);
        printf("  ");
        */
        arb_printn(I3, 20, ARB_STR_MORE);
        printf("   total ");
        arb_printn(I, 20, ARB_STR_MORE);
        printf("\n");
    }

    arb_div(I,I,pi,prec);
    arb_printn(I, 20, ARB_STR_MORE);
    acb_poly_clear(P);
    arb_clear(sigma);
    arb_clear(a);
    arb_clear(b); arb_clear(pi);
    arb_clear(I);
    arb_clear(I2);
    arb_clear(I3);
    mag_clear(tol);
    mag_clear(err);
}

