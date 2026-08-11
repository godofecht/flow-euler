// Project Euler 807: Loops of Ropes
// P(n) = A(2n-1, n) / (2n-1)!  where A(m,k) is an Eulerian number (OEIS A025585).
// a(n) = sum_{j=0..n-1} (-1)^j * C(2n, j) * (n-j)^(2n-1)
// Exact big-integer arithmetic via GMP, then round half-up to 10 decimals.
#include <gmp.h>

double p807_native(void) {
    unsigned long n = 80;
    unsigned long m = 2 * n - 1;  /* 159 */

    mpz_t total, term, binom, power;
    mpz_init(total);
    mpz_init(term);
    mpz_init(binom);
    mpz_init(power);

    for (unsigned long j = 0; j < n; j++) {
        unsigned long base = n - j;
        mpz_bin_uiui(binom, 2 * n, j);
        mpz_ui_pow_ui(power, base, m);
        mpz_mul(term, binom, power);
        if (j & 1UL) {
            mpz_sub(total, total, term);
        } else {
            mpz_add(total, total, term);
        }
    }

    mpz_t den;
    mpz_init(den);
    mpz_fac_ui(den, m);  /* (2n-1)! */

    /* Round half-up to 10 decimal places: q = (total * 10^10) / den. */
    mpz_t scaled, q, r, twice_r;
    mpz_init(scaled);
    mpz_init(q);
    mpz_init(r);
    mpz_init(twice_r);

    mpz_mul_ui(scaled, total, 10000000000UL);
    mpz_tdiv_qr(q, r, scaled, den);
    mpz_mul_ui(twice_r, r, 2);
    if (mpz_cmp(twice_r, den) >= 0) {
        mpz_add_ui(q, q, 1);
    }

    /* q holds floor(P*1e10 + 0.5) as an exact integer (~1e9, fits in double). */
    double result = mpz_get_d(q) / 1e10;

    mpz_clear(total);
    mpz_clear(term);
    mpz_clear(binom);
    mpz_clear(power);
    mpz_clear(den);
    mpz_clear(scaled);
    mpz_clear(q);
    mpz_clear(r);
    mpz_clear(twice_r);

    return result;
}
