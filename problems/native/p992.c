/* Project Euler 992
   Journey counting on a path graph modulo a prime.
   Port of /tmp/pes_ref/solvers/992.py
*/
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef long long i64;

static const i64 MOD = 987898789LL;

static i64 mod_pow(i64 a, i64 e, i64 mod) {
    i64 r = 1 % mod;
    a %= mod;
    if (a < 0) a += mod;
    while (e > 0) {
        if (e & 1) r = (i64)((__int128)r * a % mod);
        a = (i64)((__int128)a * a % mod);
        e >>= 1;
    }
    return r;
}

static i64 mod_inv(i64 a, i64 mod) { return mod_pow(a, mod - 2, mod); }

static i64 *g_fact;
static i64 *g_inv_fact;
static i64 g_limit;

static void build_combinatorics(i64 limit, i64 mod) {
    g_limit = limit;
    g_fact = malloc((size_t)(limit + 1) * sizeof(i64));
    g_inv_fact = malloc((size_t)(limit + 1) * sizeof(i64));
    g_fact[0] = 1 % mod;
    for (i64 i = 1; i <= limit; i++)
        g_fact[i] = (i64)((__int128)g_fact[i - 1] * i % mod);
    g_inv_fact[limit] = mod_inv(g_fact[limit], mod);
    for (i64 i = limit; i >= 1; i--)
        g_inv_fact[i - 1] = (i64)((__int128)g_inv_fact[i] * i % mod);
}

static i64 comb(i64 n, i64 r) {
    if (r < 0 || r > n) return 0;
    return (i64)((__int128)g_fact[n] * g_inv_fact[r] % MOD * g_inv_fact[n - r] % MOD);
}

static i64 endpoint_count(i64 n, i64 k, i64 end) {
    if (n == 0) return 1;

    /* right[i] for 0 <= i < n */
    i64 *right = malloc((size_t)n * sizeof(i64));
    right[0] = k - (end == 0 ? 1 : 0);
    if (n >= 2)
        right[1] = 2 - (end == 1 ? 1 : 0);
    for (i64 i = 2; i < n; i++)
        right[i] = 1 + right[i - 2] - (end == i ? 1 : 0);

    i64 ways = 1;
    for (i64 v = 1; v < n; v++) {
        i64 out_degree = k + v - (end == v ? 1 : 0);
        if (v < end) {
            ways = (i64)((__int128)ways * comb(out_degree - 1, right[v] - 1) % MOD);
        } else if (v == end) {
            ways = (i64)((__int128)ways * comb(out_degree, right[v]) % MOD);
        } else {
            ways = (i64)((__int128)ways * comb(out_degree - 1, right[v]) % MOD);
        }
    }
    free(right);
    return ways;
}

static i64 journey_count(i64 n, i64 k) {
    i64 total = 0;
    for (i64 end = 0; end <= n; end++) {
        total += endpoint_count(n, k, end);
        if (total >= MOD) total -= MOD;
    }
    return total;
}

long long p992_native(void) {
    i64 n = 500;
    i64 ks[5] = {1, 10, 100, 1000, 10000};
    i64 max_k = 0;
    for (int i = 0; i < 5; i++)
        if (ks[i] > max_k) max_k = ks[i];
    build_combinatorics(max_k + n, MOD);

    i64 answer = 0;
    for (int i = 0; i < 5; i++) {
        answer += journey_count(n, ks[i]);
        if (answer >= MOD) answer -= MOD;
    }

    free(g_fact);
    free(g_inv_fact);
    return answer;
}
