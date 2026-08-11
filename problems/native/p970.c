/*
 * Project Euler 970 - Kangaroo Hopping over Sixes
 *
 * H(n) = 2n + 2/3 + eps, where eps = 2*Re(e^{lam*n}/lam) is a tiny
 * correction with lam = 1 + W_1(-1/e) (Lambert W, branch 1).
 *
 * For n = 10^6 the asymptotic formula is exact to far more digits than
 * needed.  We precompute lam to 70 digits, then use MPFR at ~256-bit
 * precision for the n-dependent arithmetic (phase reduction mod 2*pi,
 * log10 of the correction, digit extraction).
 */

#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>

long long p970_native(void) {
    mpfr_set_default_prec(256);

    /* lam = 1 + W_1(-1/e), precomputed to 70 decimal digits */
    mpfr_t lam_re, lam_im, lam_abs, lam_arg;
    mpfr_init_set_str(lam_re,
        "-2.088843015613043855957086716774947500545693741036729673239112544244607",
        10, MPFR_RNDN);
    mpfr_init_set_str(lam_im,
        "7.461489285654254556906116612186415334509094993202209240934411391411877",
        10, MPFR_RNDN);
    mpfr_init_set_str(lam_abs,
        "7.748360310659838754659859240216219375966049464762612291981467582861770",
        10, MPFR_RNDN);
    mpfr_init_set_str(lam_arg,
        "1.843758551210239598129985611715645443454302414441256154232297011153285",
        10, MPFR_RNDN);

    mpfr_t n, pi, two_pi, neg_pi, theta, c, a, ln10, tmp, tmp2, s, frac;
    mpfr_inits(n, pi, two_pi, neg_pi, theta, c, a, ln10,
               tmp, tmp2, s, frac, (mpfr_ptr)0);

    mpfr_set_ui(n, 1000000, MPFR_RNDN);
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
    mpfr_neg(neg_pi, pi, MPFR_RNDN);

    /* ln(10) for converting between log10 and ln */
    mpfr_set_ui(tmp, 10, MPFR_RNDN);
    mpfr_log(ln10, tmp, MPFR_RNDN);

    /* theta = Im(lam)*n - arg(lam), reduced to [-pi, pi] */
    mpfr_mul(theta, lam_im, n, MPFR_RNDN);
    mpfr_sub(theta, theta, lam_arg, MPFR_RNDN);
    mpfr_fmod(theta, theta, two_pi, MPFR_RNDN);
    if (mpfr_cmp(theta, pi) > 0)
        mpfr_sub(theta, theta, two_pi, MPFR_RNDN);
    else if (mpfr_cmp(theta, neg_pi) < 0)
        mpfr_add(theta, theta, two_pi, MPFR_RNDN);

    /* c = cos(theta) */
    mpfr_cos(c, theta, MPFR_RNDN);

    /* a = log10(2/|lam|) + Re(lam)*n/ln(10) + log10|cos(theta)| */
    mpfr_set_ui(tmp, 2, MPFR_RNDN);
    mpfr_div(tmp, tmp, lam_abs, MPFR_RNDN);
    mpfr_log10(tmp, tmp, MPFR_RNDN);
    mpfr_set(a, tmp, MPFR_RNDN);

    mpfr_mul(tmp, lam_re, n, MPFR_RNDN);
    mpfr_div(tmp, tmp, ln10, MPFR_RNDN);
    mpfr_add(a, a, tmp, MPFR_RNDN);

    mpfr_abs(tmp, c, MPFR_RNDN);
    mpfr_log10(tmp, tmp, MPFR_RNDN);
    mpfr_add(a, a, tmp, MPFR_RNDN);

    /* L = floor(-a) */
    mpfr_neg(tmp, a, MPFR_RNDN);
    mpfr_floor(tmp, tmp);
    long long L = mpfr_get_si(tmp, MPFR_RNDN);

    /* delta = sign(c) * 10^(a + L) */
    mpfr_add_si(tmp, a, (long)L, MPFR_RNDN);   /* a + L */
    mpfr_mul(tmp, tmp, ln10, MPFR_RNDN);       /* (a+L)*ln(10) */
    mpfr_exp(tmp, tmp, MPFR_RNDN);             /* 10^(a+L) */
    if (mpfr_sgn(c) < 0)
        mpfr_neg(tmp, tmp, MPFR_RNDN);

    /* s = 2/3 + delta, frac = s - floor(s) */
    mpfr_set_ui(s, 2, MPFR_RNDN);
    mpfr_div_ui(s, s, 3, MPFR_RNDN);
    mpfr_add(s, s, tmp, MPFR_RNDN);
    mpfr_floor(frac, s);
    mpfr_sub(frac, s, frac, MPFR_RNDN);

    /* Extract first 8 digits after decimal point that are not 6 */
    long long result = 0;
    int count = 0;
    for (int i = 0; i < 200 && count < 8; i++) {
        mpfr_mul_ui(frac, frac, 10, MPFR_RNDN);
        mpfr_floor(tmp, frac);
        long long d = mpfr_get_si(tmp, MPFR_RNDN);
        mpfr_sub(frac, frac, tmp, MPFR_RNDN);
        if (d != 6) {
            result = result * 10 + d;
            count++;
        }
    }

    mpfr_clears(lam_re, lam_im, lam_abs, lam_arg, n, pi, two_pi, neg_pi,
                theta, c, a, ln10, tmp, tmp2, s, frac, (mpfr_ptr)0);
    return result;
}
