/* Project Euler 796: A Grand Shuffle
 *
 * Expected stopping time for drawing cards without replacement from 10 decks
 * (540 cards) until every suit, rank, and deck design has appeared.
 *
 * Uses inclusion-exclusion over missing suits (a), ranks (b), and deck
 * designs (c). For each combination the allowed card count is
 *   M(a,b,c) = (4-a)*(13-b)*(10-c) + 2*(10-c)
 * and the coefficient is (-1)^(a+b+c+1) * C(4,a)*C(13,b)*C(10,c).
 *
 * E[T] = sum_M coeff[M] * S(M, N)
 * where S(M,N) = sum_{k=0}^{N-1} C(M,k)/C(N,k), computed via the recurrence
 *   r_0 = 1, r_k = r_{k-1} * (M-k+1)/(N-k+1).
 *
 * MPFR is used for 80+ digit precision to match the reference Decimal solver.
 */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>

static void sum_choose_ratios(mpfr_t result, long M, long N,
                              mpfr_t r, mpfr_t tmp) {
    /* S(M,N) = sum_{k=0}^{N-1} C(M,k)/C(N,k) */
    mpfr_set_ui(result, 1, MPFR_RNDN);  /* s = 1 (k=0 term) */
    mpfr_set_ui(r, 1, MPFR_RNDN);       /* r = 1 */

    long max_k = (N - 1 < M) ? (N - 1) : M;
    for (long k = 1; k <= max_k; k++) {
        mpfr_mul_ui(r, r, (unsigned long)(M - k + 1), MPFR_RNDN);
        mpfr_div_ui(r, r, (unsigned long)(N - k + 1), MPFR_RNDN);
        mpfr_add(result, result, r, MPFR_RNDN);
    }
}

double p796_native(void) {
    mpfr_set_default_prec(280);  /* ~84 decimal digits */

    const long N = 540;
    const int max_M = 540;

    mpfr_t *coeffs = malloc((max_M + 1) * sizeof(mpfr_t));
    for (int i = 0; i <= max_M; i++) mpfr_init_set_ui(coeffs[i], 0, MPFR_RNDN);

    mpz_t coef_val;
    mpz_init(coef_val);

    for (int a = 0; a <= 4; a++) {
        for (int b = 0; b <= 13; b++) {
            for (int c = 0; c <= 10; c++) {
                if (a == 0 && b == 0 && c == 0) continue;

                /* coef = (-1)^(a+b+c+1) * C(4,a) * C(13,b) * C(10,c) */
                mpz_set_ui(coef_val, 1);
                mpz_mul_ui(coef_val, coef_val, 1);
                /* Compute C(4,a) * C(13,b) * C(10,c) */
                mpz_t tmpz;
                mpz_init(tmpz);
                mpz_bin_uiui(tmpz, 4, a);
                mpz_mul(coef_val, coef_val, tmpz);
                mpz_bin_uiui(tmpz, 13, b);
                mpz_mul(coef_val, coef_val, tmpz);
                mpz_bin_uiui(tmpz, 10, c);
                mpz_mul(coef_val, coef_val, tmpz);
                mpz_clear(tmpz);

                if ((a + b + c) % 2 == 0) mpz_neg(coef_val, coef_val);

                int decks_left = 10 - c;
                long M = (long)(4 - a) * (13 - b) * decks_left + 2 * decks_left;

                mpfr_add_z(coeffs[M], coeffs[M], coef_val, MPFR_RNDN);
            }
        }
    }
    mpz_clear(coef_val);

    mpfr_t total, s_val, r, tmp;
    mpfr_init_set_ui(total, 0, MPFR_RNDN);
    mpfr_init(s_val);
    mpfr_init(r);
    mpfr_init(tmp);

    for (int M = 0; M <= max_M; M++) {
        if (mpfr_zero_p(coeffs[M])) continue;
        sum_choose_ratios(s_val, M, N, r, tmp);
        mpfr_mul(tmp, coeffs[M], s_val, MPFR_RNDN);
        mpfr_add(total, total, tmp, MPFR_RNDN);
    }

    /* Round to 8 decimal places (ROUND_HALF_UP) */
    mpfr_t rounded;
    mpfr_init(rounded);
    mpfr_mul_ui(rounded, total, 100000000, MPFR_RNDN);
    mpfr_round(rounded, rounded);
    mpfr_div_ui(rounded, rounded, 100000000, MPFR_RNDN);

    double result = mpfr_get_d(rounded, MPFR_RNDN);

    for (int i = 0; i <= max_M; i++) mpfr_clear(coeffs[i]);
    free(coeffs);
    mpfr_clear(total); mpfr_clear(s_val); mpfr_clear(r);
    mpfr_clear(tmp); mpfr_clear(rounded);

    return result;
}
