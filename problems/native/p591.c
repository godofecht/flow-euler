
#include <mpfr.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define PREC 500  /* bits ~ 150 decimal digits */

static void frac_mpfr(mpfr_t out, const mpfr_t x) {
    mpfr_t fl; mpfr_init2(fl, PREC);
    mpfr_floor(fl, x);
    mpfr_sub(out, x, fl, MPFR_RNDN);
    mpfr_clear(fl);
}

static long ceil_mpfr_long(const mpfr_t x) {
    mpfr_t c; mpfr_init2(c, PREC);
    mpfr_ceil(c, x);
    long v = mpfr_get_si(c, MPFR_RNDN);
    mpfr_clear(c);
    return v;
}

static long nearest_mpfr_long(const mpfr_t x) {
    mpfr_t r; mpfr_init2(r, PREC);
    mpfr_round(r, x);
    long v = mpfr_get_si(r, MPFR_RNDN);
    mpfr_clear(r);
    return v;
}

static void chudnovsky_pi(mpfr_t pi) {
    /* Use MPFR's const_pi for correctness/speed; algorithmically π is well-defined. */
    mpfr_const_pi(pi, MPFR_RNDN);
}

static int is_square(int n) {
    int r = (int)floor(sqrt((double)n) + 0.5);
    return r * r == n;
}

static int sqrt_cf_period(int D, int *a0_out, int *period, int *plen) {
    int a0 = (int)floor(sqrt((double)D));
    *a0_out = a0;
    long m = 0, d = 1, a = a0;
    int n = 0;
    for (;;) {
        m = d * a - m;
        d = (D - m * m) / d;
        a = (a0 + m) / d;
        period[n++] = (int)a;
        if (a == 2 * a0) break;
        if (n > 10000) return -1;
    }
    *plen = n;
    return 0;
}

static long best_b_positive(mpfr_t alpha, mpfr_t beta, long B, const int *period, int plen) {
    if (B <= 0) return 0;
    /* dynamic arrays */
    int cap = 64;
    int *a = calloc((size_t)cap, sizeof(int));
    long *q = calloc((size_t)cap, sizeof(long));
    mpfr_t *delta = calloc((size_t)cap, sizeof(mpfr_t));
    for (int i = 0; i < cap; i++) mpfr_init2(delta[i], PREC);

    a[0] = 0;
    q[0] = 1;
    long q_minus1 = 0;
    mpfr_set(delta[0], alpha, MPFR_RNDN);
    mpfr_t delta_minus1; mpfr_init2(delta_minus1, PREC); mpfr_set_ui(delta_minus1, 1, MPFR_RNDN);

    int k = 1;
    int extra = 6;
    int len = 1; /* number of filled indices including 0 */
    while (1) {
        if (k + 2 >= cap) {
            int ncap = cap * 2;
            a = realloc(a, (size_t)ncap * sizeof(int));
            q = realloc(q, (size_t)ncap * sizeof(long));
            delta = realloc(delta, (size_t)ncap * sizeof(mpfr_t));
            for (int i = cap; i < ncap; i++) mpfr_init2(delta[i], PREC);
            memset(a + cap, 0, (size_t)(ncap - cap) * sizeof(int));
            memset(q + cap, 0, (size_t)(ncap - cap) * sizeof(long));
            cap = ncap;
        }
        int ak = period[(k - 1) % plen];
        a[k] = ak;
        long qk = (long)ak * q[k - 1] + q_minus1;
        q_minus1 = q[k - 1];
        q[k] = qk;

        mpfr_t tmp; mpfr_init2(tmp, PREC);
        if (k == 1) {
            mpfr_mul_si(tmp, delta[0], -ak, MPFR_RNDN);
            mpfr_add(delta[k], tmp, delta_minus1, MPFR_RNDN);
        } else {
            mpfr_mul_si(tmp, delta[k - 1], -ak, MPFR_RNDN);
            mpfr_add(delta[k], tmp, delta[k - 2], MPFR_RNDN);
        }
        mpfr_clear(tmp);
        len = k + 1;
        if (qk > B) {
            extra--;
            if (extra <= 0) break;
        }
        k++;
        if (k > 500) break;
    }
    int max_i = len - 1;

    int *b_digits = calloc((size_t)(max_i + 1), sizeof(int));
    mpfr_t beta_rem; mpfr_init2(beta_rem, PREC);
    mpfr_set(beta_rem, beta, MPFR_RNDN);
    for (int i = 1; i <= max_i; i++) {
        mpfr_t ratio; mpfr_init2(ratio, PREC);
        mpfr_div(ratio, beta_rem, delta[i - 1], MPFR_RNDN);
        long bi = ceil_mpfr_long(ratio);
        if (bi > a[i]) bi = a[i];
        if (bi < 0) bi = 0;
        b_digits[i] = (int)bi;
        mpfr_mul_si(ratio, delta[i - 1], bi, MPFR_RNDN);
        mpfr_sub(beta_rem, ratio, beta_rem, MPFR_RNDN);
        mpfr_clear(ratio);
    }

    long *prefix = calloc((size_t)(max_i + 1), sizeof(long));
    long s = 0;
    for (int i = 1; i <= max_i; i++) {
        s += (long)b_digits[i] * q[i - 1];
        prefix[i] = s;
    }

    /* candidates */
    long *cand_r = malloc(sizeof(long) * 100000);
    long *cand_l = malloc(sizeof(long) * 100000);
    int nr = 0, nl = 0;
    cand_r[nr++] = 0;
    cand_l[nl++] = 0;
    for (int k2 = 1; 2 * k2 < max_i + 1; k2++) {
        int idx_even = 2 * k2;
        int idx_odd = 2 * k2 - 1;
        if (idx_even > max_i) break;
        long P = prefix[idx_odd];
        long step = q[idx_odd];
        for (int j = 0; j < b_digits[idx_even]; j++) {
            long n = P + (long)j * step;
            if (0 <= n && n <= B) cand_r[nr++] = n;
        }
    }
    for (int k2 = 0; 2 * k2 + 1 <= max_i; k2++) {
        int idx = 2 * k2;
        int idx_next = idx + 1;
        if (idx_next > max_i) break;
        long P = prefix[idx];
        long step = q[idx];
        for (int j = 0; j < b_digits[idx_next]; j++) {
            long n = P + (long)j * step;
            if (0 <= n && n <= B) cand_l[nl++] = n;
        }
    }

    mpfr_t best_r_gap, best_l_gap, x, gap, one;
    mpfr_inits2(PREC, best_r_gap, best_l_gap, x, gap, one, (mpfr_ptr)0);
    mpfr_set_ui(best_r_gap, 1, MPFR_RNDN);
    mpfr_set_ui(best_l_gap, 1, MPFR_RNDN);
    mpfr_set_ui(one, 1, MPFR_RNDN);
    long best_r_n = 0, best_l_n = 0;
    for (int i = 0; i < nr; i++) {
        long n = cand_r[i];
        mpfr_mul_si(x, alpha, n, MPFR_RNDN);
        frac_mpfr(x, x);
        mpfr_sub(gap, x, beta, MPFR_RNDN);
        if (mpfr_sgn(gap) < 0) mpfr_add(gap, gap, one, MPFR_RNDN);
        if (mpfr_cmp(gap, best_r_gap) < 0) { mpfr_set(best_r_gap, gap, MPFR_RNDN); best_r_n = n; }
    }
    for (int i = 0; i < nl; i++) {
        long n = cand_l[i];
        mpfr_mul_si(x, alpha, n, MPFR_RNDN);
        frac_mpfr(x, x);
        mpfr_sub(gap, beta, x, MPFR_RNDN);
        if (mpfr_sgn(gap) < 0) mpfr_add(gap, gap, one, MPFR_RNDN);
        if (mpfr_cmp(gap, best_l_gap) < 0) { mpfr_set(best_l_gap, gap, MPFR_RNDN); best_l_n = n; }
    }
    long ans = (mpfr_cmp(best_r_gap, best_l_gap) < 0) ? best_r_n : best_l_n;

    for (int i = 0; i < cap; i++) mpfr_clear(delta[i]);
    free(delta); free(a); free(q); free(b_digits); free(prefix); free(cand_r); free(cand_l);
    mpfr_clear(delta_minus1); mpfr_clears(best_r_gap, best_l_gap, x, gap, one, beta_rem, (mpfr_ptr)0);
    return ans;
}

static void bqa_pi_d(int d, long n, mpfr_t pi, mpfr_t beta_pi, long *a_out, long *b_out) {
    mpfr_t sqrt_d, alpha; mpfr_inits2(PREC, sqrt_d, alpha, (mpfr_ptr)0);
    mpfr_set_ui(sqrt_d, (unsigned)d, MPFR_RNDN);
    mpfr_sqrt(sqrt_d, sqrt_d, MPFR_RNDN);
    int a0, period[4096], plen;
    sqrt_cf_period(d, &a0, period, &plen);
    mpfr_sub_ui(alpha, sqrt_d, (unsigned)a0, MPFR_RNDN);

    mpfr_t tmp; mpfr_init2(tmp, PREC);
    mpfr_add_si(tmp, pi, n, MPFR_RNDN);
    mpfr_div(tmp, tmp, sqrt_d, MPFR_RNDN);
    mpfr_floor(tmp, tmp);
    long Bpos = mpfr_get_si(tmp, MPFR_RNDN);
    if (Bpos < 0) Bpos = 0; if (Bpos > n) Bpos = n;

    mpfr_si_sub(tmp, n, pi, MPFR_RNDN);
    mpfr_div(tmp, tmp, sqrt_d, MPFR_RNDN);
    mpfr_floor(tmp, tmp);
    long Bneg = mpfr_get_si(tmp, MPFR_RNDN);
    if (Bneg < 0) Bneg = 0; if (Bneg > n) Bneg = n;

    long b_pos = best_b_positive(alpha, beta_pi, Bpos, period, plen);
    mpfr_t one_m_beta; mpfr_init2(one_m_beta, PREC);
    mpfr_ui_sub(one_m_beta, 1, beta_pi, MPFR_RNDN);
    long t = best_b_positive(alpha, one_m_beta, Bneg, period, plen);
    long b_neg = -t;

    mpfr_t a_real, err1, err2, ab; mpfr_inits2(PREC, a_real, err1, err2, ab, (mpfr_ptr)0);
    /* candidate b_pos */
    mpfr_mul_si(ab, sqrt_d, b_pos, MPFR_RNDN);
    mpfr_sub(a_real, pi, ab, MPFR_RNDN);
    long a1 = nearest_mpfr_long(a_real);
    if (a1 > n) a1 = n; if (a1 < -n) a1 = -n;
    mpfr_mul_si(tmp, sqrt_d, b_pos, MPFR_RNDN);
    mpfr_add_si(ab, tmp, a1, MPFR_RNDN);
    mpfr_sub(err1, ab, pi, MPFR_RNDN);
    mpfr_abs(err1, err1, MPFR_RNDN);

    mpfr_mul_si(tmp, sqrt_d, b_neg, MPFR_RNDN);
    mpfr_sub(a_real, pi, tmp, MPFR_RNDN);
    long a2 = nearest_mpfr_long(a_real);
    if (a2 > n) a2 = n; if (a2 < -n) a2 = -n;
    mpfr_mul_si(tmp, sqrt_d, b_neg, MPFR_RNDN);
    mpfr_add_si(ab, tmp, a2, MPFR_RNDN);
    mpfr_sub(err2, ab, pi, MPFR_RNDN);
    mpfr_abs(err2, err2, MPFR_RNDN);

    if (mpfr_cmp(err2, err1) < 0) { *a_out = a2; *b_out = b_neg; }
    else { *a_out = a1; *b_out = b_pos; }

    mpfr_clears(sqrt_d, alpha, tmp, one_m_beta, a_real, err1, err2, ab, (mpfr_ptr)0);
}

long long pe591_answer(void) {
    mpfr_set_default_prec(PREC);
    mpfr_t pi, beta; mpfr_inits2(PREC, pi, beta, (mpfr_ptr)0);
    chudnovsky_pi(pi);
    frac_mpfr(beta, pi);
    long long total = 0;
    long n = 10000000000000LL; /* 10^13 */
    for (int d = 2; d < 100; d++) {
        if (is_square(d)) continue;
        long a, b;
        bqa_pi_d(d, n, pi, beta, &a, &b);
        total += (a >= 0 ? a : -a);
    }
    mpfr_clears(pi, beta, (mpfr_ptr)0);
    return total;
}
