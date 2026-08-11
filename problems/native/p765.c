/* Project Euler 765 - Trillionaire
 *
 * Under a fair-coin measure the wealth process is a martingale, so the number
 * of paths on which we can guarantee terminal wealth >= M is at most
 * floor(2^n / M).  Spend that budget on the most-likely paths (most wins).
 *
 * Uses GMP for exact big-integer arithmetic.
 */

#include <gmp.h>
#include <stdlib.h>
#include <string.h>

double p765_native(void) {
    const int n = 1000;
    /* M = 10^12 */
    mpz_t M;
    mpz_init_set_ui(M, 10UL);
    mpz_pow_ui(M, M, 12);

    /* total_paths = 2^n */
    mpz_t total_paths;
    mpz_init(total_paths);
    mpz_set_ui(total_paths, 1UL);
    mpz_mul_2exp(total_paths, total_paths, n);

    /* budget_paths = total_paths // M */
    mpz_t budget_paths;
    mpz_init(budget_paths);
    mpz_fdiv_q(budget_paths, total_paths, M);

    /* binomial coefficients C(n,k) for k=0..n */
    mpz_t *comb = malloc((n + 1) * sizeof(mpz_t));
    for (int k = 0; k <= n; k++) mpz_init(comb[k]);
    mpz_set_ui(comb[0], 1UL);
    for (int k = 1; k <= n; k++) {
        /* C(n,k) = C(n,k-1) * (n-k+1) / k */
        mpz_mul_ui(comb[k], comb[k - 1], (unsigned long)(n - k + 1));
        mpz_fdiv_q_ui(comb[k], comb[k], (unsigned long)k);
    }

    /* suffix[k] = number of paths with >= k wins = sum_{j>=k} C(n,j) */
    mpz_t *suffix = malloc((n + 2) * sizeof(mpz_t));
    for (int k = 0; k <= n + 1; k++) mpz_init(suffix[k]);
    mpz_set_ui(suffix[n + 1], 0UL);
    for (int k = n; k >= 0; k--) {
        mpz_add(suffix[k], suffix[k + 1], comb[k]);
    }

    /* Find k0 such that suffix[k0] > budget_paths >= suffix[k0+1] */
    int k0 = -1;
    for (int k = n; k >= 0; k--) {
        if (mpz_cmp(suffix[k], budget_paths) > 0 &&
            mpz_cmp(budget_paths, suffix[k + 1]) >= 0) {
            k0 = k;
            break;
        }
    }

    /* rem = budget_paths - suffix[k0+1] */
    mpz_t rem;
    mpz_init(rem);
    mpz_sub(rem, budget_paths, suffix[k0 + 1]);

    /* Precompute pow2[i] = 2^i, pow3[i] = 3^i */
    mpz_t *pow2 = malloc((n + 1) * sizeof(mpz_t));
    mpz_t *pow3 = malloc((n + 1) * sizeof(mpz_t));
    for (int i = 0; i <= n; i++) { mpz_init(pow2[i]); mpz_init(pow3[i]); }
    mpz_set_ui(pow2[0], 1UL);
    mpz_set_ui(pow3[0], 1UL);
    for (int i = 1; i <= n; i++) {
        mpz_mul_2exp(pow2[i], pow2[i - 1], 1);
        mpz_mul_ui(pow3[i], pow3[i - 1], 3UL);
    }

    /* den = 5^n */
    mpz_t den;
    mpz_init(den);
    mpz_set_ui(den, 5UL);
    mpz_pow_ui(den, den, n);

    /* num = sum over k>=k0+1 of C(n,k)*3^k*2^(n-k)  +  rem*3^k0*2^(n-k0) */
    mpz_t num, term;
    mpz_init(num);
    mpz_init(term);
    mpz_set_ui(num, 0UL);

    for (int k = k0 + 1; k <= n; k++) {
        /* term = comb[k] * pow3[k] * pow2[n-k] */
        mpz_mul(term, comb[k], pow3[k]);
        mpz_mul(term, term, pow2[n - k]);
        mpz_add(num, num, term);
    }

    if (mpz_sgn(rem) > 0) {
        mpz_mul(term, rem, pow3[k0]);
        mpz_mul(term, term, pow2[n - k0]);
        mpz_add(num, num, term);
    }

    /* Round num/den to 10 decimal places (round half up). */
    const int digits = 10;
    mpz_t scale;
    mpz_init(scale);
    mpz_set_ui(scale, 10UL);
    mpz_pow_ui(scale, scale, digits);

    /* q = (num * scale) / den, r = remainder */
    mpz_t scaled_num, q, r;
    mpz_init(scaled_num);
    mpz_init(q);
    mpz_init(r);
    mpz_mul(scaled_num, num, scale);
    mpz_fdiv_qr(q, r, scaled_num, den);

    /* round half up: if 2*r >= den, q += 1 */
    mpz_mul_2exp(r, r, 1);
    if (mpz_cmp(r, den) >= 0) {
        mpz_add_ui(q, q, 1UL);
    }

    /* q now holds the value * 10^10.  Convert to double. */
    /* result = q / 10^10 as a double */
    double result = mpz_get_d(q) / 1e10;

    /* cleanup */
    mpz_clear(M); mpz_clear(total_paths); mpz_clear(budget_paths);
    for (int k = 0; k <= n; k++) mpz_clear(comb[k]);
    free(comb);
    for (int k = 0; k <= n + 1; k++) mpz_clear(suffix[k]);
    free(suffix);
    mpz_clear(rem);
    for (int i = 0; i <= n; i++) { mpz_clear(pow2[i]); mpz_clear(pow3[i]); }
    free(pow2); free(pow3);
    mpz_clear(den); mpz_clear(num); mpz_clear(term);
    mpz_clear(scale); mpz_clear(scaled_num); mpz_clear(q); mpz_clear(r);

    return result;
}
