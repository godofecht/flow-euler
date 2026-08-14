// Project Euler 833: Square Triangle Products
// S(10^35) mod 136101521.
// Uses GMP for exact big integer arithmetic.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#define MOD 136101521LL

// GCD
static long long gcd_ll(long long a, long long b) {
    while (b) { long long t = a % b; a = b; b = t; }
    return a;
}

// Lucas U-sequence: U_0=0, U_1=1, U_k = P*U_{k-1} - U_{k-2}
// Returns (U_i, U_j) as mpz_t values
static void lucas_U_pair_mpz(long long n, int i, int j, mpz_t ui, mpz_t uj) {
    long long P = 4 * n + 2;
    int kmax = i > j ? i : j;
    mpz_t u0, u1, tmp;
    mpz_init_set_ui(u0, 0);
    mpz_init_set_ui(u1, 1);
    mpz_init(tmp);

    if (kmax == 0) {
        mpz_set_ui(ui, 0);
        mpz_set_ui(uj, 0);
        mpz_clears(u0, u1, tmp, NULL);
        return;
    }
    if (kmax == 1) {
        mpz_set_ui(ui, i == 1 ? 1 : 0);
        mpz_set_ui(uj, j == 1 ? 1 : 0);
        mpz_clears(u0, u1, tmp, NULL);
        return;
    }

    if (i == 1) mpz_set(ui, u1); else mpz_set_ui(ui, 0);
    if (j == 1) mpz_set(uj, u1); else mpz_set_ui(uj, 0);

    for (int k = 2; k <= kmax; k++) {
        // u1 = P*u1 - u0
        mpz_mul_ui(tmp, u1, P);
        mpz_sub(tmp, tmp, u0);
        mpz_set(u0, u1);
        mpz_set(u1, tmp);
        if (k == i) mpz_set(ui, u1);
        if (k == j) mpz_set(uj, u1);
    }

    mpz_clears(u0, u1, tmp, NULL);
}

// c_value(n, i, j) = T(n) * U_i * U_j where T(n) = n*(n+1)/2
static void c_value_mpz(long long n, int i, int j, mpz_t result) {
    mpz_t ui, uj, t;
    mpz_init(ui);
    mpz_init(uj);
    mpz_init(t);

    lucas_U_pair_mpz(n, i, j, ui, uj);

    // T(n) = n*(n+1)/2
    mpz_set_ui(t, n);
    mpz_mul_ui(t, t, n + 1);
    mpz_divexact_ui(t, t, 2);

    mpz_mul(result, t, ui);
    mpz_mul(result, result, uj);

    mpz_clears(ui, uj, t, NULL);
}

// Compare c_value(n, i, j) with N (as mpz_t). Returns -1, 0, 1.
static int c_value_cmp(long long n, int i, int j, const mpz_t N) {
    mpz_t cv;
    mpz_init(cv);
    c_value_mpz(n, i, j, cv);
    int cmp = mpz_cmp(cv, N);
    mpz_clear(cv);
    return cmp;
}

// Find largest n such that c_value(n, i, j) <= N
static long long max_n_for_pair(int i, int j, const mpz_t N) {
    long long lo = 0, hi = 1;
    while (c_value_cmp(hi, i, j, N) <= 0)
        hi *= 2;
    lo = hi / 2;

    while (lo + 1 < hi) {
        long long mid = (lo + hi) / 2;
        if (c_value_cmp(mid, i, j, N) <= 0)
            lo = mid;
        else
            hi = mid;
    }
    return lo;
}

// Forward differences: given f(0..d), return Newton coefficients a_m
static void forward_differences(mpz_t *values, int d, mpz_t *coeffs) {
    mpz_t *cur = malloc((d + 1) * sizeof(mpz_t));
    for (int k = 0; k <= d; k++) {
        mpz_init_set(cur[k], values[k]);
    }

    int n = d + 1;
    int ci = 0;
    while (n > 0) {
        mpz_init_set(coeffs[ci], cur[0]);
        ci++;
        for (int k = 0; k < n - 1; k++) {
            mpz_sub(cur[k], cur[k + 1], cur[k]);
        }
        n--;
    }

    for (int k = 0; k <= d; k++) mpz_clear(cur[k]);
    free(cur);
}

// Exact binomial C(n, k)
static void binom_exact(long long n, int k, mpz_t result) {
    if (k < 0 || k > n) { mpz_set_ui(result, 0); return; }
    if (k == 0 || k == n) { mpz_set_ui(result, 1); return; }
    if (k > n - k) k = n - k;
    mpz_set_ui(result, 1);
    for (int t = 1; t <= k; t++) {
        mpz_mul_ui(result, result, n - k + t);
        mpz_divexact_ui(result, result, t);
    }
}

// sum_{n=0..M} f(n) using Newton series
static void sum_poly_prefix(mpz_t *fvals, int deg, long long M, mpz_t result) {
    mpz_t *a = malloc((deg + 1) * sizeof(mpz_t));
    for (int i = 0; i <= deg; i++) mpz_init(a[i]);

    forward_differences(fvals, deg, a);

    mpz_t binom, term;
    mpz_init(binom);
    mpz_init(term);
    mpz_set_ui(result, 0);

    long long Mp1 = M + 1;
    for (int m = 0; m <= deg; m++) {
        binom_exact(Mp1, m + 1, binom);
        mpz_mul(term, a[m], binom);
        mpz_add(result, result, term);
    }

    for (int i = 0; i <= deg; i++) mpz_clear(a[i]);
    free(a);
    mpz_clears(binom, term, NULL);
}

// Determine max j such that c(1,1,j) <= N
// For n=1, T(1)=1, P=6, c(1,1,j) = U_j(6)
static int maxJ_for_N(const mpz_t N) {
    long long P = 6;
    mpz_t u0, u1, tmp;
    mpz_init_set_ui(u0, 0);
    mpz_init_set_ui(u1, 1);
    mpz_init(tmp);

    int maxj = 1;
    for (int k = 2; k < 400; k++) {
        mpz_mul_ui(tmp, u1, P);
        mpz_sub(tmp, tmp, u0);
        mpz_set(u0, u1);
        mpz_set(u1, tmp);
        if (mpz_cmp(u1, N) > 0) break;
        maxj = k;
    }

    mpz_clears(u0, u1, tmp, NULL);
    return maxj;
}

long long p833_native(void) {
    mpz_t N;
    mpz_init(N);
    mpz_ui_pow_ui(N, 10, 35);  // N = 10^35

    int maxj = maxJ_for_N(N);

    mpz_t total, cv, fval, sum_result;
    mpz_init_set_ui(total, 0);
    mpz_init(cv);
    mpz_init(fval);
    mpz_init(sum_result);

    for (int i = 1; i < maxj; i++) {
        for (int j = i + 1; j <= maxj; j++) {
            if (gcd_ll(i, j) != 1) continue;

            // Check c_value(1, i, j) > N
            c_value_mpz(1, i, j, cv);
            if (mpz_cmp(cv, N) > 0) continue;

            long long M = max_n_for_pair(i, j, N);
            if (M <= 0) continue;

            int deg = i + j;
            mpz_t *fvals = malloc((deg + 1) * sizeof(mpz_t));
            for (int n = 0; n <= deg; n++) {
                mpz_init(fvals[n]);
                c_value_mpz(n, i, j, fvals[n]);
            }

            sum_poly_prefix(fvals, deg, M, sum_result);
            mpz_add(total, total, sum_result);

            for (int n = 0; n <= deg; n++) mpz_clear(fvals[n]);
            free(fvals);
        }
    }

    long long result = mpz_fdiv_ui(total, MOD);

    mpz_clears(N, total, cv, fval, sum_result, NULL);
    return result;
}
