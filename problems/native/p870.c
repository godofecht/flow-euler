#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

/* Project Euler 870: Stone Game IV
   Transition values T(i) computed via g-base generation with exact integer arithmetic.
   T(1) = 1, and each subsequent transition is the minimum upper endpoint b_k/b_{i-1}.
   Target: T(123456) as a double.

   The g-base values grow exponentially, so we use GMP for exact arithmetic. */

#define MAX_STEPS 3000

/* Compute the next transition value after q = q_num/q_den.
   Returns the result as a reduced fraction via out_num/out_den (mpz). */
static void next_transition(const mpz_t q_num, const mpz_t q_den,
                            mpz_t out_num, mpz_t out_den) {
    static mpz_t b[MAX_STEPS + 2];
    static int initialized = 0;
    if (!initialized) {
        for (int i = 0; i < MAX_STEPS + 2; i++) mpz_init(b[i]);
        initialized = 1;
    }

    int len = 0;
    mpz_set_ui(b[len], 1);
    len++;
    int i_ptr = 0;

    mpz_t best_num, best_den, rhs, tmp_cmp1, tmp_cmp2, bk_plus_bi;
    mpz_init(best_num); mpz_init(best_den);
    mpz_init(rhs); mpz_init(tmp_cmp1); mpz_init(tmp_cmp2);
    mpz_init(bk_plus_bi);

    mpz_set_ui(best_num, 0);
    mpz_set_ui(best_den, 1);
    int have_best = 0;

    for (int step = 0; step < MAX_STEPS; step++) {
        /* bk = b[len-1] */
        mpz_mul(rhs, b[len - 1], q_den);  /* rhs = bk * q_den */

        /* advance i_ptr until q_num*b[i_ptr] >= rhs */
        while (i_ptr < len) {
            mpz_mul(tmp_cmp1, q_num, b[i_ptr]);
            if (mpz_cmp(tmp_cmp1, rhs) >= 0) break;
            i_ptr++;
        }
        if (i_ptr >= len)
            break;

        /* candidate upper endpoint: b_k / b[i_ptr - 1] */
        if (i_ptr > 0) {
            if (!have_best) {
                mpz_set(best_num, b[len - 1]);
                mpz_set(best_den, b[i_ptr - 1]);
                have_best = 1;
            } else {
                /* compare b_k * best_den < best_num * b[i_ptr-1] */
                mpz_mul(tmp_cmp1, b[len - 1], best_den);
                mpz_mul(tmp_cmp2, best_num, b[i_ptr - 1]);
                if (mpz_cmp(tmp_cmp1, tmp_cmp2) < 0) {
                    mpz_set(best_num, b[len - 1]);
                    mpz_set(best_den, b[i_ptr - 1]);
                }
            }
        }

        /* extend: b[len] = bk + b[i_ptr] */
        mpz_add(bk_plus_bi, b[len - 1], b[i_ptr]);
        mpz_set(b[len], bk_plus_bi);
        len++;
    }

    if (!have_best) {
        mpz_set_ui(out_num, 1);
        mpz_set_ui(out_den, 1);
    } else {
        mpz_t g;
        mpz_init(g);
        mpz_gcd(g, best_num, best_den);
        mpz_divexact(out_num, best_num, g);
        mpz_divexact(out_den, best_den, g);
        mpz_clear(g);
    }

    mpz_clear(best_num); mpz_clear(best_den);
    mpz_clear(rhs); mpz_clear(tmp_cmp1); mpz_clear(tmp_cmp2);
    mpz_clear(bk_plus_bi);
}

double p870_native(void) {
    mpz_t q_num, q_den, nn, dd;
    mpz_init(q_num); mpz_init(q_den);
    mpz_init(nn); mpz_init(dd);

    mpz_set_ui(q_num, 1);
    mpz_set_ui(q_den, 1); /* T(1) = 1 */

    for (int i = 1; i < 123456; i++) {
        next_transition(q_num, q_den, nn, dd);
        mpz_set(q_num, nn);
        mpz_set(q_den, dd);
    }

    /* Convert to double */
    double result = mpz_get_d(q_num) / mpz_get_d(q_den);

    mpz_clear(q_num); mpz_clear(q_den);
    mpz_clear(nn); mpz_clear(dd);

    return result;
}
