// Project Euler 881: Divisor Graph Width
// Find smallest n with g(n) >= 10000, where g(n) is the maximum coefficient
// of the product of (1 + x + ... + x^{e_i}) for the exponent vector of n.
// Branch-and-bound DFS over non-increasing exponent sequences.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int PRIMES[400];
static int num_primes;

static void sieve_primes(int limit) {
    char *sieve = calloc(limit + 1, 1);
    memset(sieve, 1, limit + 1);
    if (limit >= 0) { sieve[0] = 0; sieve[1] = 0; }
    for (int p = 2; (long long)p * p <= limit; p++) {
        if (sieve[p]) {
            for (int m = p * p; m <= limit; m += p)
                sieve[m] = 0;
        }
    }
    num_primes = 0;
    for (int i = 2; i <= limit; i++)
        if (sieve[i]) PRIMES[num_primes++] = i;
    free(sieve);
}

// Convolve poly with (1 + x + ... + x^e), sliding window sum.
static int *convolve_with_ones(const int *poly, int poly_len, int e, int *out_len) {
    *out_len = poly_len + e;
    int *out = calloc(*out_len, sizeof(int));
    long long *pref = calloc(poly_len + 1, sizeof(long long));
    long long s = 0;
    for (int i = 0; i < poly_len; i++) {
        s += poly[i];
        pref[i + 1] = s;
    }
    for (int j = 0; j < *out_len; j++) {
        int lo = j - e;
        if (lo < 0) lo = 0;
        int hi = j;
        if (hi >= poly_len) hi = poly_len - 1;
        if (lo <= hi)
            out[j] = (int)(pref[hi + 1] - pref[lo]);
    }
    free(pref);
    return out;
}

static int max_array(const int *arr, int len) {
    int m = arr[0];
    for (int i = 1; i < len; i++)
        if (arr[i] > m) m = arr[i];
    return m;
}

static __int128 best_n;

static void dfs(int idx, int prev_e, const int *poly, int poly_len, int peak,
                __int128 n, int target) {
    if (n >= best_n) return;
    if (peak >= target) {
        best_n = n;
        return;
    }
    if (idx >= num_primes) return;

    int p = PRIMES[idx];

    // Maximum exponent allowed by monotonicity and keeping n < best_n
    int max_e = 0;
    __int128 t = n;
    while (max_e < prev_e) {
        t *= p;
        if (t >= best_n) break;
        max_e++;
    }
    if (max_e == 0) return;

    // Precompute powers
    __int128 powers[64];
    powers[0] = 1;
    for (int e = 1; e <= max_e; e++)
        powers[e] = powers[e-1] * p;

    for (int e = max_e; e >= 1; e--) {
        __int128 n2 = n * powers[e];
        if (n2 >= best_n) continue;
        int out_len;
        int *poly2 = convolve_with_ones(poly, poly_len, e, &out_len);
        int peak2 = max_array(poly2, out_len);
        dfs(idx + 1, e, poly2, out_len, peak2, n2, target);
        free(poly2);
    }
}

static __int128 initial_upper_bound(int target) {
    int *poly = malloc(sizeof(int));
    poly[0] = 1;
    int poly_len = 1;
    __int128 n = 1;
    int i = 0;
    while (1) {
        int out_len;
        int *poly2 = convolve_with_ones(poly, poly_len, 1, &out_len);
        free(poly);
        poly = poly2;
        poly_len = out_len;
        n *= PRIMES[i];
        i++;
        if (max_array(poly, poly_len) >= target) {
            free(poly);
            return n;
        }
    }
}

long long p881_native(void) {
    sieve_primes(400);

    int target = 10000;
    best_n = initial_upper_bound(target);

    int *initial_poly = malloc(sizeof(int));
    initial_poly[0] = 1;
    dfs(0, 60, initial_poly, 1, 1, 1, target);
    free(initial_poly);

    return (long long)best_n;
}
