/*
    Copyright (C) 2018 Association des collaborateurs de D.H.J Polymath

    This is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_calc.h"

static void
_acb_get_rad_ubound_mag(mag_t res, const acb_t z, slong prec)
{
    arf_t err;
    arf_init(err);
    acb_get_rad_ubound_arf(err, z, prec);
    arf_get_mag(res, err);
    arf_clear(err);
}

typedef struct
{
    acb_t z;
    arb_t t;
    arb_t T;
} thirdbound_integrand_param_struct;
typedef thirdbound_integrand_param_struct thirdbound_integrand_param_t[1];
typedef thirdbound_integrand_param_struct *thirdbound_integrand_param_ptr;
typedef const thirdbound_integrand_param_struct *thirdbound_integrand_param_srcptr;

static void
thirdbound_integrand_param_init(thirdbound_integrand_param_t p,
        const acb_t z, const arb_t t, const arb_t T)
{
    acb_init(p->z);
    arb_init(p->t);
    arb_init(p->T);
    acb_set(p->z, z);
    arb_set(p->t, t);
    arb_set(p->T, T);
}

static void
thirdbound_integrand_param_clear(thirdbound_integrand_param_t p)
{
    acb_clear(p->z);
    arb_clear(p->t);
    arb_clear(p->T);
}

void thirdbound_Q(arb_t res, const arb_t t, const arb_t y,
        const arb_t X, const arb_t T, slong prec)
{
    arb_t z1, z2, z3;

    arb_init(z1);
    arb_mul_2exp_si(z1, X, 1);
    arb_add(z1, z1, y, prec);
    arb_sub_ui(z1, z1, 1, prec);
    arb_div(z1, z1, t, prec);

    arb_init(z2);
    arb_set_ui(z2, 3);
    arb_div(z2, z2, T, prec);
    arb_mul_2exp_si(z2, z2, -2);

    arb_init(z3);
    arb_log_hypot(z3, X, T, prec);
    arb_mul_2exp_si(z3, z3, -1);

    arb_const_log_sqrt2pi(res, prec);
    arb_add(res, res, z1, prec);
    arb_sub(res, res, z2, prec);
    arb_sub(res, res, z3, prec);

    arb_clear(z1);
    arb_clear(z2);
    arb_clear(z3);
}

void thirdbound_F(acb_t res, const acb_t s, slong prec)
{
    acb_t sd2, sm1;
    acb_t z1, z2, z3;

    acb_init(sd2);
    acb_mul_2exp_si(sd2, s, -1);

    acb_init(sm1);
    acb_sub_ui(sm1, s, 1, prec);

    acb_init(z1);
    acb_log(z1, sm1, prec);

    acb_init(z2);
    acb_const_pi(z2, prec);
    acb_log(z2, z2, prec);
    acb_mul(z2, z2, sd2, prec);

    acb_init(z3);
    acb_log(z3, sd2, prec);
    acb_mul(z3, z3, sm1, prec);
    acb_mul_2exp_si(z3, z3, -1);

    acb_log(res, s, prec);
    acb_add(res, res, z1, prec);
    acb_sub(res, res, z2, prec);
    acb_add(res, res, z3, prec);
    acb_sub(res, res, sd2, prec);

    acb_clear(sd2);
    acb_clear(sm1);

    acb_clear(z1);
    acb_clear(z2);
    acb_clear(z3);
}

void riemann_xi(acb_t res, const acb_t s, slong prec)
{
    acb_t pi, z1, z2, z3, z4;

    acb_init(pi);
    acb_const_pi(pi, prec);

    acb_init(z1);
    acb_sub_ui(z1, s, 1, prec);
    acb_mul(z1, z1, s, prec);
    acb_mul_2exp_si(z1, z1, -1);

    acb_init(z2);
    acb_mul_2exp_si(z2, s, -1);
    acb_neg(z2, z2);
    acb_pow(z2, pi, z2, prec);

    acb_init(z3);
    acb_mul_2exp_si(z3, s, -1);
    acb_gamma(z3, z3, prec);

    acb_init(z4);
    acb_zeta(z4, s, prec);

    acb_mul(res, z1, z2, prec);
    acb_mul(res, res, z3, prec);
    acb_mul(res, res, z4, prec);

    acb_clear(pi);
    acb_clear(z1);
    acb_clear(z2);
    acb_clear(z3);
    acb_clear(z4);
}




static void
_thirdbound_remainder_helper(arb_t res,
        const arb_t t, const acb_t z,
        const arb_t X, const arb_t T, slong prec)
{
    arb_t z1, z2, z3, z4;
    acb_t u;
    arb_srcptr x, y;

    x = acb_realref(z);
    y = acb_imagref(z);

    arb_init(z1);
    arb_mul_2exp_si(z1, x, -1);
    arb_sub(z1, T, z1, prec);
    arb_sqr(z1, z1, prec);
    arb_div(z1, z1, t, prec);
    arb_mul_2exp_si(z1, z1, -2);

    arb_init(z2);
    arb_set_str(z2, "0.33", prec);
    arb_sub(z2, T, z2, prec);
    arb_mul_ui(z2, z2, 6, prec);
    arb_inv(z2, z2, prec);

    acb_init(u);
    acb_set_arb_arb(u, X, T);
    thirdbound_F(u, u, prec);

    arb_init(z3);
    arb_set(z3, acb_realref(u));

    arb_init(z4);
    arb_sub_ui(z4, y, 1, prec);
    arb_mul_2exp_si(z4, z4, -1);
    arb_add(z4, X, z4, prec);
    arb_sqr(z4, z4, prec);
    arb_div(z4, z4, t, prec);

    arb_add(res, z1, z2, prec);
    arb_add(res, res, z3, prec);
    arb_sub(res, res, z4, prec);
    arb_exp(res, res, prec);

    arb_clear(z1);
    arb_clear(z2);
    arb_clear(z3);
    arb_clear(z4);

    acb_clear(u);
}


void thirdbound_H_remainder(arb_t res,
        const arb_t t, const acb_t z,
        const arb_t X, const arb_t T, slong prec)
{
    arb_t q, u;
    arb_srcptr x, y;

    x = acb_realref(z);
    y = acb_imagref(z);

    arb_init(q);
    thirdbound_Q(q, t, y, X, T, prec);

    arb_init(u);
    arb_mul_2exp_si(u, t, 5);
    arb_sqrt(u, u, prec);
    arb_mul(u, u, q, prec);

    _thirdbound_remainder_helper(res, t, z, X, T, prec);

    arb_div(res, res, u, prec);

    arb_clear(q);
    arb_clear(u);
}


void thirdbound_Hprime_remainder(arb_t res,
        const arb_t t, const acb_t z,
        const arb_t X, const arb_t T, slong prec)
{
    arb_t qinv, u;
    arb_t z1, z2, w;
    arb_srcptr x, y;

    x = acb_realref(z);
    y = acb_imagref(z);

    arb_init(qinv);
    thirdbound_Q(qinv, t, y, X, T, prec);
    arb_inv(qinv, qinv, prec);

    arb_init(u);
    arb_pow_ui(u, t, 3, prec);
    arb_mul_2exp_si(u, u, 5);
    arb_sqrt(u, u, prec);

    arb_init(z1);
    arb_mul_2exp_si(z1, x, -1);
    arb_sub(z1, T, z1, prec);
    arb_abs(z1, z1);

    arb_init(z2);
    arb_sub_ui(z2, y, 1, prec);
    arb_mul_2exp_si(z2, z2, -1);
    arb_add(z2, z2, X, prec);
    
    arb_init(w);
    arb_add(w, z1, z2, prec);
    arb_add(w, w, qinv, prec);
    arb_mul(w, w, qinv, prec);

    _thirdbound_remainder_helper(res, t, z, X, T, prec);
    arb_mul(res, res, w, prec);
    arb_div(res, res, u, prec);

    arb_clear(qinv);
    arb_clear(u);
    arb_clear(z1);
    arb_clear(z2);
    arb_clear(w);
}




static void
_thirdbound_integrand_helper(acb_t res1, acb_t res2,
        const acb_t s, const acb_t z, const arb_t T, slong prec)
{
    acb_t onei, iT;
    acb_t z1, z2;
    arb_srcptr x, y;

    x = acb_realref(z);
    y = acb_imagref(z);

    acb_init(onei);
    acb_onei(onei);

    acb_init(iT);
    acb_mul_arb(iT, onei, T, prec);

    acb_init(z1);
    acb_add(z1, s, iT, prec);
    riemann_xi(z1, z1, prec);

    acb_init(z2);
    acb_mul_arb(z2, onei, x, prec);
    acb_sub_arb(z2, z2, y, prec);
    acb_add_ui(z2, z2, 1, prec);
    acb_mul_2exp_si(z2, z2, -1);
    acb_sub(z2, iT, z2, prec);

    acb_swap(res1, z1);
    acb_swap(res2, z2);

    acb_clear(z1);
    acb_clear(z2);

    acb_clear(onei);
    acb_clear(iT);
}



void thirdbound_H_integrand(acb_t res, const acb_t z,
        const arb_t t, const acb_t s, const arb_t T, slong prec)
{
    acb_t h1, h2, z1, z2, w1, w2;

    acb_init(h1);
    acb_init(h2);
    _thirdbound_integrand_helper(h1, h2, s, z, T, prec);

    acb_init(z1);
    acb_add(z1, s, h2, prec);
    acb_sqr(z1, z1, prec);
    acb_neg(z1, z1);
    acb_div_arb(z1, z1, t, prec);
    acb_exp(z1, z1, prec);

    acb_init(z2);
    acb_sub(z2, h2, s, prec);
    acb_add_ui(z2, z2, 1, prec);
    acb_sqr(z2, z2, prec);
    acb_neg(z2, z2);
    acb_div_arb(z2, z2, t, prec);
    acb_exp(z2, z2, prec);

    acb_init(w1);
    acb_mul(w1, h1, z1, prec);

    acb_init(w2);
    acb_conj(w2, h1);
    acb_mul(w2, w2, z2, prec);

    acb_add(res, w1, w2, prec);

    acb_clear(h1);
    acb_clear(h2);
    acb_clear(z1);
    acb_clear(z2);
    acb_clear(w1);
    acb_clear(w2);
}

void thirdbound_Hprime_integrand(acb_t res, const acb_t z,
        const arb_t t, const acb_t s, const arb_t T, slong prec)
{
    acb_t h1, h2, u1, u2, z1, z2, w1, w2;

    acb_init(h1);
    acb_init(h2);
    _thirdbound_integrand_helper(h1, h2, s, z, T, prec);

    acb_init(u1);
    acb_add(u1, s, h2, prec);

    acb_init(u2);
    acb_add_ui(u2, h2, 1, prec);
    acb_sub(u2, u2, s, prec);

    acb_init(z1);
    acb_sqr(z1, u1, prec);
    acb_neg(z1, z1);
    acb_div_arb(z1, z1, t, prec);
    acb_exp(z1, z1, prec);

    acb_init(z2);
    acb_sqr(z2, u2, prec);
    acb_neg(z2, z2);
    acb_div_arb(z2, z2, t, prec);
    acb_exp(z2, z2, prec);

    acb_init(w1);
    acb_mul(w1, h1, u1, prec);
    acb_mul(w1, w1, z1, prec);

    acb_init(w2);
    acb_conj(w2, h1);
    acb_mul(w2, w2, u2, prec);
    acb_mul(w2, w2, z2, prec);

    acb_add(res, w1, w2, prec);

    acb_clear(h1);
    acb_clear(h2);
    acb_clear(u1);
    acb_clear(u2);
    acb_clear(z1);
    acb_clear(z2);
    acb_clear(w1);
    acb_clear(w2);
}

static int
f_thirdbound_H_integrand(acb_ptr res, const acb_t u, void * param,
        slong order, slong prec)
{
    thirdbound_integrand_param_srcptr p;

    if (order > 1)
        flint_abort();

    p = param;
    thirdbound_H_integrand(res, p->z, p->t, u, p->T, prec);

    return 0;
}

static int
f_thirdbound_Hprime_integrand(acb_ptr res, const acb_t u, void * param,
        slong order, slong prec)
{
    thirdbound_integrand_param_srcptr p;

    if (order > 1)
        flint_abort();

    p = param;
    thirdbound_Hprime_integrand(res, p->z, p->t, u, p->T, prec);

    return 0;
}

void
thirdbound_H_proper_integral(acb_t res,
        const mag_t abs_estimate,
        const arb_t lower_limit, const arb_t upper_limit,
        const acb_t z, const arb_t t, const arb_t T, slong prec)
{
    thirdbound_integrand_param_t p;
    acb_t a_lim, b_lim;
    acb_calc_integrate_opt_t options;
    int calc_result;
    mag_t abs_tol;
    slong rel_goal = prec;

    thirdbound_integrand_param_init(p, z, t, T);

    acb_init(a_lim);
    acb_init(b_lim);
    acb_set_arb(a_lim, lower_limit);
    acb_set_arb(b_lim, upper_limit);

    mag_init(abs_tol);
    mag_mul_2exp_si(abs_tol, abs_estimate, -prec);

    acb_calc_integrate_opt_init(options);
    options->verbose = 0; 
    calc_result = acb_calc_integrate(
            res, f_thirdbound_H_integrand, p, a_lim, b_lim,
            rel_goal, abs_tol, options, prec);
    if (calc_result == ARB_CALC_NO_CONVERGENCE)
    {
        acb_indeterminate(res);
    }
    else
    {
        arb_t u;
        arb_init(u);
        arb_const_pi(u, prec);
        arb_mul(u, u, t, prec);
        arb_sqrt(u, u, prec);
        acb_div_arb(res, res, u, prec);
        acb_mul_2exp_si(res, res, -3);
        arb_clear(u);
    }

    thirdbound_integrand_param_clear(p);
    acb_clear(a_lim);
    acb_clear(b_lim);
    mag_clear(abs_tol);
}


void
thirdbound_Hprime_proper_integral(acb_t res,
        const mag_t abs_estimate,
        const arb_t lower_limit, const arb_t upper_limit,
        const acb_t z, const arb_t t, const arb_t T, slong prec)
{
    thirdbound_integrand_param_t p;
    acb_t a_lim, b_lim;
    acb_calc_integrate_opt_t options;
    int calc_result;
    mag_t abs_tol;
    slong rel_goal = prec;

    thirdbound_integrand_param_init(p, z, t, T);

    acb_init(a_lim);
    acb_init(b_lim);
    acb_set_arb(a_lim, lower_limit);
    acb_set_arb(b_lim, upper_limit);

    mag_init(abs_tol);
    mag_mul_2exp_si(abs_tol, abs_estimate, -prec);

    acb_calc_integrate_opt_init(options);
    options->verbose = 0; 
    calc_result = acb_calc_integrate(
            res, f_thirdbound_Hprime_integrand, p, a_lim, b_lim,
            rel_goal, abs_tol, options, prec);
    if (calc_result == ARB_CALC_NO_CONVERGENCE)
    {
        acb_indeterminate(res);
    }
    else
    {
        fmpq_t q;
        arb_t sqrtpi, tpow;
        acb_t onei;

        fmpq_init(q);
        fmpq_set_si(q, 3, 2);

        arb_init(sqrtpi);
        arb_const_sqrt_pi(sqrtpi, prec);

        arb_init(tpow);
        arb_pow_fmpq(tpow, t, q, prec);

        acb_init(onei);
        acb_onei(onei);

        acb_mul(res, res, onei, prec);
        acb_div_arb(res, res, sqrtpi, prec);
        acb_div_arb(res, res, tpow, prec);
        acb_mul_2exp_si(res, res, -3);

        fmpq_clear(q);
        arb_clear(sqrtpi);
        arb_clear(tpow);
        acb_clear(onei);
    }

    thirdbound_integrand_param_clear(p);
    acb_clear(a_lim);
    acb_clear(b_lim);
    mag_clear(abs_tol);
}


int thirdbound_H_remainder_is_permissible(const arb_t t, const arb_t y,
        const arb_t X, const arb_t T, slong prec)
{
    arb_t q, zero, one, two, half;
    int result;

    arb_init(q);
    arb_init(zero);
    arb_init(one);
    arb_init(two);
    arb_init(half);

    arb_zero(zero);
    arb_one(one);
    arb_set_ui(two, 2);
    arb_mul_2exp_si(half, one, -1);
    thirdbound_Q(q, t, y, X, T, prec);

    result = (
            arb_gt(X, one) &&
            arb_ge(y, zero) &&
            arb_ge(T, two) &&
            arb_le(t, half) &&
            arb_gt(q, zero));

    arb_clear(q);
    arb_clear(zero);
    arb_clear(one);
    arb_clear(two);
    arb_clear(half);

    return result;
}

void
thirdbound_H_integral(acb_t res,
        const mag_t m_estimate,
        const acb_t z, const arb_t t, const arb_t T, slong prec)
{
    arb_t lower, upper, remainder;
    acb_t value, partial;
    mag_t m_partial_sum, tail_mag_lb, tail_mag_ub;
    slong r, p, q, tmp;

    arb_init(lower);
    arb_init(upper);
    arb_init(remainder);

    acb_init(value);
    acb_init(partial);

    arb_init(remainder);

    mag_init(m_partial_sum);
    mag_init(tail_mag_lb);
    mag_init(tail_mag_ub);

    r = 0;
    p = 0;
    q = 1;
    while (1)
    {
        tmp = p + q;
        arb_set_ui(lower, 1 + r);
        arb_set_ui(upper, 1 + r + q);
        arb_mul_2exp_si(lower, lower, -1);
        arb_mul_2exp_si(upper, upper, -1);
        
        r = r + q;
        p = q;
        q = tmp;

        thirdbound_H_proper_integral(value,
                m_estimate, lower, upper, z, t, T, prec);
        acb_add(partial, partial, value, prec);

        if (thirdbound_H_remainder_is_permissible(
                    t, acb_imagref(z), upper, T, prec))
        {
            thirdbound_H_remainder(remainder,
                    t, z, upper, T, prec);
            arb_get_mag(tail_mag_ub, remainder);

            _acb_get_rad_ubound_mag(m_partial_sum, partial, prec);
            arb_get_mag_lower(tail_mag_lb, remainder);

            if (mag_cmp(tail_mag_lb, m_partial_sum) <= 0)
            {
                acb_set(res, partial);
                acb_add_error_mag(res, tail_mag_ub);
                break;
            }
        }
    }

    arb_clear(lower);
    arb_clear(upper);

    acb_clear(value);
    acb_clear(partial);

    arb_clear(remainder);

    mag_clear(m_partial_sum);
    mag_clear(tail_mag_lb);
    mag_clear(tail_mag_ub);
}


void
thirdbound_Hprime_integral(acb_t res,
        const mag_t m_estimate,
        const acb_t z, const arb_t t, const arb_t T, slong prec)
{
    arb_t lower, upper, remainder;
    acb_t value, partial;
    mag_t m_partial_sum, tail_mag_lb, tail_mag_ub;
    slong r, p, q, tmp;

    arb_init(lower);
    arb_init(upper);
    arb_init(remainder);

    acb_init(value);
    acb_init(partial);

    arb_init(remainder);

    mag_init(m_partial_sum);
    mag_init(tail_mag_lb);
    mag_init(tail_mag_ub);

    r = 0;
    p = 0;
    q = 1;
    while (1)
    {
        tmp = p + q;
        arb_set_ui(lower, 1 + r);
        arb_set_ui(upper, 1 + r + q);
        arb_mul_2exp_si(lower, lower, -1);
        arb_mul_2exp_si(upper, upper, -1);
        
        r = r + q;
        p = q;
        q = tmp;

        thirdbound_Hprime_proper_integral(value,
                m_estimate, lower, upper, z, t, T, prec);
        acb_add(partial, partial, value, prec);

        if (thirdbound_H_remainder_is_permissible(
                    t, acb_imagref(z), upper, T, prec))
        {
            thirdbound_Hprime_remainder(remainder,
                    t, z, upper, T, prec);
            arb_get_mag(tail_mag_ub, remainder);

            _acb_get_rad_ubound_mag(m_partial_sum, partial, prec);
            arb_get_mag_lower(tail_mag_lb, remainder);

            if (mag_cmp(tail_mag_lb, m_partial_sum) <= 0)
            {
                acb_set(res, partial);
                acb_add_error_mag(res, tail_mag_ub);
                break;
            }
        }
    }

    arb_clear(lower);
    arb_clear(upper);

    acb_clear(value);
    acb_clear(partial);

    arb_clear(remainder);

    mag_clear(m_partial_sum);
    mag_clear(tail_mag_lb);
    mag_clear(tail_mag_ub);
}

static void
H_approximation(arb_t res, const arb_t x, const arb_t y, slong prec)
{
    if (arb_contains_nonpositive(x))
    {
        arb_one(res);
    }
    else
    {
        arb_t logx, pixd8;

        arb_init(logx);
        arb_log(logx, x, prec);

        arb_init(pixd8);
        arb_const_pi(pixd8, prec);
        arb_mul(pixd8, pixd8, x, prec);
        arb_mul_2exp_si(pixd8, pixd8, -3);

        arb_add_si(res, y, 7, prec);
        arb_mul(res, res, logx, prec);
        arb_mul_2exp_si(res, res, -2);
        arb_sub(res, res, pixd8, prec);
        arb_exp(res, res, prec);

        arb_clear(logx);
        arb_clear(pixd8);
    }
}

int main(int argc, char *argv[])
{
    int result;
    arb_t t, T;
    acb_t z, h;
    arb_t half, four;
    slong digits, prec, high_prec;
    arb_t estimate;
    mag_t m_estimate;
    slong n;
    arb_ptr x, y;

    arb_init(t);
    arb_init(T);

    acb_init(z);
    acb_init(h);

    arb_init(half);
    arb_init(four);

    arb_init(estimate);
    mag_init(m_estimate);

    x = acb_realref(z);
    y = acb_imagref(z);

    arb_one(half);
    arb_mul_2exp_si(half, half, -1);
    arb_set_ui(four, 4);
    
    result = EXIT_SUCCESS;

    if (argc != 6)
    {
        result = EXIT_FAILURE;
        goto finish;
    }

    n = atol(argv[4]);
    if (n != 0 && n != 1)
    {
        result = EXIT_FAILURE;
        goto finish;
    }

    digits = atol(argv[5]);
    high_prec = digits * 3.32192809488736 + 10;

    prec = 16;
    acb_indeterminate(h);
    while (arb_rel_accuracy_bits(acb_realref(h)) < high_prec ||
           arb_rel_accuracy_bits(acb_imagref(h)) < high_prec)
    {
        prec *= 2;
        arb_set_str(t, argv[1], prec);
        arb_set_str(acb_realref(z), argv[2], prec);
        arb_set_str(acb_imagref(z), argv[3], prec);

        if (!arb_is_nonnegative(x) ||
            !arb_is_nonnegative(y) ||
            !arb_le(t, half))
        {
            result = EXIT_FAILURE;
            goto finish;
        }

        arb_const_pi(T, prec);
        arb_mul_2exp_si(T, T, -1);
        arb_add(T, T, x, prec);
        arb_mul_2exp_si(T, T, -1);
        arb_max(T, T, four, prec);

        H_approximation(estimate, x, y, prec);
        arb_get_mag(m_estimate, estimate);

        if (n == 0)
        {
            thirdbound_H_integral(h, m_estimate, z, t, T, prec);
            if (arb_is_zero(x) || arb_is_zero(y))
            {
                arb_zero(acb_imagref(h));
            }
        }
        else if (n == 1)
        {
            thirdbound_Hprime_integral(h, m_estimate, z, t, T, prec);
            if (arb_is_zero(x))
            {
                arb_zero(acb_realref(h));
            }
            if (arb_is_zero(y))
            {
                arb_zero(acb_imagref(h));
            }
        }
        else
        {
            flint_abort();
        }

    }

    flint_printf("Re: ");
    arf_printd(arb_midref(acb_realref(h)), digits);
    flint_printf("\n");
    flint_printf("Im: ");
    arf_printd(arb_midref(acb_imagref(h)), digits);
    flint_printf("\n");

finish:

    if (result == EXIT_FAILURE)
    {
        flint_printf("Usage:\n");
        flint_printf("%s t x y n d\n\n", argv[0]);
        flint_printf(
                "Evaluate the nth derivative of H_t(x + yi) "
                "to d significant digits "
                "where H is the function involved "
                "in the definition of the De Bruijn-Newman constant.\n"
                "Requires x >= 0 and y >= 0 and t <= 1/2 and n in {0, 1}.\n");
    }

    arb_clear(t);
    arb_clear(T);
    arb_clear(half);
    arb_clear(four);

    acb_clear(z);
    acb_clear(h);

    arb_clear(estimate);
    mag_clear(m_estimate);

    flint_cleanup();
    return result;
}
