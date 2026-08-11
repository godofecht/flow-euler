// Project Euler 708: Twos Are All You Need
// S(N) = sum_{n<=N} 2^{Omega(n)}, computed for N = 10^14.
//
// Identity: 2^{Omega(n)} = sum_{d|n, d powerful} g(d), where a powerful number
// d = prod p_i^{e_i} (e_i >= 2) carries weight g(d) = prod 2^{e_i - 2}.
// Then S(N) = sum_{d powerful, d<=N} g(d) * D(N/d), where
// D(x) = sum_{k<=x} tau(k) = #{(a,b): a*b <= x} = 2*sum_{i<=s} floor(x/i) - s^2,
// s = isqrt(x). Small x use a precomputed prefix of tau; large x are cached.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef long long i64;
typedef unsigned long long u64;

static i64 isqrt_i64(i64 n) {
    if (n <= 0) return 0;
    i64 x = (i64)sqrt((double)n);
    while (x > 0 && x * x > n) x--;
    while ((x + 1) * (x + 1) <= n) x++;
    return x;
}

// Odd-only sieve; returns malloc'd prime list (2, 3, 5, ...).
static int *primes_upto(int limit, int *count_out) {
    if (limit < 2) { *count_out = 0; return NULL; }
    int size = limit / 2 + 1;
    char *sieve = malloc(size);
    memset(sieve, 1, size);
    sieve[0] = 0; // 1 is not prime
    int r = (int)isqrt_i64(limit);
    for (int p = 3; p <= r; p += 2) {
        if (sieve[p / 2]) {
            int start = p * p;
            for (int j = start; j <= limit; j += 2 * p) {
                sieve[j / 2] = 0;
            }
        }
    }
    int cnt = 1; // 2
    for (int i = 1; i < size; i++) if (sieve[i]) cnt++;
    int *primes = malloc((size_t)cnt * sizeof(int));
    primes[0] = 2;
    int k = 1;
    for (int i = 1; i < size; i++) if (sieve[i]) primes[k++] = 2 * i + 1;
    free(sieve);
    *count_out = cnt;
    return primes;
}

// Prefix sums of tau(k) (number of divisors) for k = 1..limit.
static i64 *build_divisor_prefix(int limit) {
    int *spf = calloc((size_t)limit + 1, sizeof(int));
    for (int i = 2; i <= limit; i++) {
        if (spf[i] == 0) {
            spf[i] = i;
            if ((i64)i * i <= limit) {
                for (int j = i * i; j <= limit; j += i) {
                    if (spf[j] == 0) spf[j] = i;
                }
            }
        }
    }
    int *d = calloc((size_t)limit + 1, sizeof(int));
    int *exp = calloc((size_t)limit + 1, sizeof(int));
    d[1] = 1;
    for (int i = 2; i <= limit; i++) {
        int p = spf[i];
        int m = i / p;
        if (m % p == 0) {
            exp[i] = exp[m] + 1;
            d[i] = d[m] / (exp[m] + 1) * (exp[i] + 1);
        } else {
            exp[i] = 1;
            d[i] = d[m] * 2;
        }
    }
    i64 *prefix = calloc((size_t)limit + 1, sizeof(i64));
    i64 total = 0;
    for (int i = 1; i <= limit; i++) {
        total += d[i];
        prefix[i] = total;
    }
    free(spf);
    free(d);
    free(exp);
    return prefix;
}

// Cache for D(x) with x > small_limit. Open-addressing hash table.
// Only ~18k distinct large values arise (powerful d < N/small_limit), so a
// small table is plenty.
#define CACHE_BITS 17
#define CACHE_SIZE (1 << CACHE_BITS)
#define CACHE_MASK (CACHE_SIZE - 1)
static i64 cache_key[CACHE_SIZE];
static i64 cache_val[CACHE_SIZE];
static char cache_used[CACHE_SIZE];

static i64 D_func(i64 x, const i64 *prefix, int small_limit) {
    if (x <= small_limit) return prefix[x];
    u64 h = (u64)x * 11400714819323198485ULL;
    int idx = (int)(h >> (64 - CACHE_BITS));
    while (cache_used[idx]) {
        if (cache_key[idx] == x) return cache_val[idx];
        idx = (idx + 1) & CACHE_MASK;
    }
    i64 s = isqrt_i64(x);
    i64 total = 0;
    i64 i = 1;
    while (i <= s) {
        i64 q = x / i;
        i64 j = x / q;
        if (j > s) j = s;
        total += q * (j - i + 1);
        i = j + 1;
    }
    i64 res = 2 * total - s * s;
    cache_used[idx] = 1;
    cache_key[idx] = x;
    cache_val[idx] = res;
    return res;
}

static i64 g_total;
static int *g_primes;
static i64 *g_prime_sq;
static int g_nprimes;
static i64 *g_prefix;
static int g_small_limit;
static i64 g_N;

static void dfs(int start_idx, i64 current_n, i64 current_g) {
    g_total += current_g * D_func(g_N / current_n, g_prefix, g_small_limit);
    i64 limit = g_N / current_n;
    for (int i = start_idx; i < g_nprimes; i++) {
        i64 p2 = g_prime_sq[i];
        if (p2 > limit) break;
        int p = g_primes[i];
        i64 pow_p = p2;
        i64 g = current_g;
        while (pow_p <= limit) {
            dfs(i + 1, current_n * pow_p, g);
            if (pow_p > limit / p) break;
            pow_p *= p;
            g <<= 1;
        }
    }
}

i64 p708_native(void) {
    i64 N = 100000000000000LL; // 10^14
    int prime_limit = (int)isqrt_i64(N); // 10^7
    int nprimes;
    int *primes = primes_upto(prime_limit, &nprimes);
    i64 *prime_sq = malloc((size_t)nprimes * sizeof(i64));
    for (int i = 0; i < nprimes; i++) prime_sq[i] = (i64)primes[i] * primes[i];

    int small_limit = 1000000;
    i64 *prefix = build_divisor_prefix(small_limit);

    g_primes = primes;
    g_prime_sq = prime_sq;
    g_nprimes = nprimes;
    g_prefix = prefix;
    g_small_limit = small_limit;
    g_N = N;
    g_total = 0;

    dfs(0, 1, 1);

    free(primes);
    free(prime_sq);
    free(prefix);
    return g_total;
}
