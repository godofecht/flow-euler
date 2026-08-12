/* Project Euler 850 - Fractions of Powers.
 * Compute floor(S(33557799775533)) mod 977676779.
 *
 * For odd k: 2*f_k(n) = n - n/m(n,k) where m(n,k) = prod p^{ceil(e/k)}.
 * 2*S(N) = sum_{odd k<=N} (sum_{n=1..N} n - H_k) = odd_count * sum_n - total_h
 * where H_k = sum_n h_k(n) and h_k(n) = n - n/m(n,k).
 *
 * H_k is computed via DFS over squarefree d with prime powers, using
 * inclusion-exclusion with Euler's totient.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef int64_t i64;
typedef __int128 i128;

#define MOD 977676779LL
#define MOD2 (2 * MOD)

static i64 g_N;
static i64 g_mod2;

static int *g_primes;
static int g_num_primes;

static struct { i64 p; int exponent; } g_factors[20];
static int g_factors_count;

static i64 g_H_ge3;

static i64 odd_count_from(i64 first, i64 limit) {
    if (first > limit) return 0;
    if (first % 2 == 0) first += 1;
    if (first > limit) return 0;
    return (limit - first) / 2 + 1;
}

static i64 imax_i64(i64 a, i64 b) { return a > b ? a : b; }
static i64 imin_i64(i64 a, i64 b) { return a < b ? a : b; }

static i64 profile_contribution(i64 d, i64 rad, i64 phi_mod, int max_exp) {
    i64 tail_first = imax_i64(3, (i64)(max_exp + 1));
    if (tail_first % 2 == 0) tail_first += 1;

    i64 total = (i64)((i128)(odd_count_from(tail_first, g_N) % g_mod2) * (g_N / (d * rad)) % g_mod2);

    i64 small_last = imin_i64(g_N, tail_first - 2);
    i64 limit = g_N / d;
    for (i64 k = 3; k <= small_last; k += 2) {
        i64 rho = 1;
        i64 denom = k - 1;
        for (int i = 0; i < g_factors_count; i++) {
            i64 p = g_factors[i].p;
            int exponent = g_factors[i].exponent;
            int e = (exponent + (int)denom - 1) / (int)denom;
            i64 p_power = 1;
            int overflow = 0;
            for (int j = 0; j < e; j++) {
                if (p_power > limit / p) { overflow = 1; break; }
                p_power *= p;
            }
            if (overflow) { rho = limit + 1; break; }
            if (rho > limit / p_power) { rho = limit + 1; break; }
            rho *= p_power;
        }
        if (rho <= limit)
            total += g_N / (d * rho);
    }
    return (i64)((i128)phi_mod * (total % g_mod2) % g_mod2);
}

static void visit(int start_idx, i64 d, i64 rad, i64 phi_mod, int max_exp) {
    g_H_ge3 = (g_H_ge3 + profile_contribution(d, rad, phi_mod, max_exp)) % g_mod2;

    for (int i = start_idx; i < g_num_primes; i++) {
        i64 p = g_primes[i];
        /* Check d * rad * p * p > N without overflow */
        if (d * rad > g_N / (p * p)) break;

        i64 next_rad = rad * p;
        i64 next_d = d * p;
        i64 phi_factor = p - 1;
        int exponent = 1;
        while (next_d <= g_N / next_rad) {
            g_factors[g_factors_count].p = p;
            g_factors[g_factors_count].exponent = exponent;
            g_factors_count++;
            visit(i + 1, next_d, next_rad,
                  (i64)((i128)phi_mod * (phi_factor % g_mod2) % g_mod2),
                  max_exp > exponent ? max_exp : exponent);
            g_factors_count--;

            exponent++;
            next_d *= p;
            phi_factor *= p;
        }
    }
}

static void sieve_primes(int limit, int **out_primes, int *out_count) {
    char *sieve = (char *)calloc((size_t)limit + 1, 1);
    for (int p = 2; (i64)p * p <= limit; p++)
        if (!sieve[p])
            for (int j = p * p; j <= limit; j += p)
                sieve[j] = 1;
    int count = 0;
    for (int p = 2; p <= limit; p++)
        if (!sieve[p]) count++;
    *out_primes = (int *)malloc((size_t)count * sizeof(int));
    int idx = 0;
    for (int p = 2; p <= limit; p++)
        if (!sieve[p]) (*out_primes)[idx++] = p;
    *out_count = count;
    free(sieve);
}

static i64 isqrt_i64(i64 n) {
    if (n < 0) return 0;
    i64 r = (i64)sqrt((double)n);
    while (r > 0 && r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r;
}

long long p850_native(void) {
    g_N = 33557799775533LL;
    g_mod2 = MOD2;
    g_factors_count = 0;
    g_H_ge3 = 0;

    i64 sqrt_N = isqrt_i64(g_N);
    sieve_primes((int)sqrt_N, &g_primes, &g_num_primes);

    visit(0, 1, 1, 1, 0);

    i64 H1 = g_N % g_mod2;
    i64 odd_count = (g_N + 1) / 2;
    i64 total_h = (H1 + g_H_ge3) % g_mod2;
    i64 sum_n = (i64)(((i128)g_N * (g_N + 1) / 2) % g_mod2);
    i64 twoS = (i64)((i128)(odd_count % g_mod2) * sum_n % g_mod2) - total_h;
    twoS %= g_mod2;
    if (twoS < 0) twoS += g_mod2;

    free(g_primes);
    return (long long)((twoS / 2) % MOD);
}
