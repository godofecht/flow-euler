/* Project Euler 929: Odd-Run Compositions
 *
 * Compute F(10^5) mod 1111124111, where F(n) counts compositions of n
 * whose maximal runs of equal parts all have odd length.
 *
 * Algorithm:
 * 1. Fibonacci numbers mod MOD up to N
 * 2. Sieve: s[n] = sum_{d|n} (-1)^(d-1) * Fib[d] mod MOD
 * 3. Invert power series f(x) = 1 - S(x) via Newton iteration
 * 4. Answer is t[N] where t = 1/f
 *
 * Uses NTT (Number Theoretic Transform) with 3 primes and CRT
 * for exact polynomial multiplication.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;
typedef __int128 i128;

#define MOD 1111124111LL
#define N 100000

/* ---- NTT primes ---- */
/* p1 = 998244353, g=3, max_ntt=2^23 */
/* p2 = 1004535809, g=3, max_ntt=2^21 */
/* p3 = 754974721, g=11, max_ntt=2^24 */
/* p1*p2*p3 ~ 7.6e26 > max_conv ~ 1.2e23 */

#define NTT_NPRIMES 3
static const i64 NTT_P[NTT_NPRIMES] = {998244353LL, 1004535809LL, 754974721LL};
static const i64 NTT_G[NTT_NPRIMES] = {3LL, 3LL, 11LL};

/* CRT precomputed: Garner's algorithm coefficients */
static i64 p1_mod_p2_inv;    /* inv(p1, p2) */
static i64 p12_mod_p3_inv;   /* inv(p1*p2 mod p3, p3) */
static i64 p1_mod_MOD;       /* p1 % MOD */
static i64 p12_mod_MOD;      /* (p1*p2) % MOD */

/* ---- Modular arithmetic ---- */

static i64 mulmod(i64 a, i64 b, i64 mod) {
    return (i64)((i128)a * b % mod);
}

static i64 powmod(i64 a, i64 e, i64 mod) {
    i64 r = 1 % mod;
    a %= mod;
    if (a < 0) a += mod;
    while (e > 0) {
        if (e & 1) r = mulmod(r, a, mod);
        a = mulmod(a, a, mod);
        e >>= 1;
    }
    return r;
}

static i64 modinv(i64 a, i64 mod) {
    return powmod(a, mod - 2, mod);
}

/* ---- NTT ---- */

static int *ntt_rev = NULL;
static int ntt_rev_n = 0;

static void ntt_bitrev(int n) {
    if (ntt_rev_n != n) {
        ntt_rev = (int *)realloc(ntt_rev, n * sizeof(int));
        ntt_rev_n = n;
    }
    int logn = 0, tmp = n;
    while (tmp > 1) { logn++; tmp >>= 1; }
    ntt_rev[0] = 0;
    for (int i = 1; i < n; i++)
        ntt_rev[i] = (ntt_rev[i >> 1] >> 1) | ((i & 1) << (logn - 1));
}

static void ntt(i64 *a, int n, int invert, i64 mod, i64 gen) {
    for (int i = 0; i < n; i++) {
        int j = ntt_rev[i];
        if (i < j) { i64 t = a[i]; a[i] = a[j]; a[j] = t; }
    }

    for (int len = 2; len <= n; len <<= 1) {
        int half = len >> 1;
        i64 w = invert ? powmod(gen, (mod - 1) - (mod - 1) / len, mod)
                       : powmod(gen, (mod - 1) / len, mod);
        for (int i = 0; i < n; i += len) {
            i64 wk = 1;
            for (int j = 0; j < half; j++) {
                i64 u = a[i + j];
                i64 v = mulmod(a[i + j + half], wk, mod);
                a[i + j] = (u + v) % mod;
                a[i + j + half] = ((u - v) % mod + mod) % mod;
                wk = mulmod(wk, w, mod);
            }
        }
    }

    if (invert) {
        i64 inv_n = modinv(n, mod);
        for (int i = 0; i < n; i++)
            a[i] = mulmod(a[i], inv_n, mod);
    }
}

/* ---- Convolution mod MOD using 3-prime NTT + CRT ---- */

static i64 *conv_mod(const i64 *a, int na, const i64 *b, int nb, int limit) {
    if (na == 0 || nb == 0) return NULL;

    int full_len = na + nb - 1;
    if (limit > full_len) limit = full_len;
    if (limit <= 0) return NULL;

    /* Small case: naive */
    if ((i64)na * nb <= 16384) {
        i64 *res = (i64 *)calloc(limit, sizeof(i64));
        for (int i = 0; i < na; i++) {
            if (a[i] == 0) continue;
            int maxj = nb;
            if (limit - i < maxj) maxj = limit - i;
            for (int j = 0; j < maxj; j++) {
                res[i + j] = (res[i + j] + a[i] * b[j]) % MOD;
                if (res[i + j] < 0) res[i + j] += MOD;
            }
        }
        return res;
    }

    int n_ntt = 1;
    while (n_ntt < full_len) n_ntt <<= 1;
    ntt_bitrev(n_ntt);

    i64 *res = (i64 *)calloc(limit, sizeof(i64));

    /* Compute NTT for each prime */
    i64 *ntt_results[NTT_NPRIMES];
    for (int pi = 0; pi < NTT_NPRIMES; pi++) {
        i64 mod = NTT_P[pi];
        i64 gen = NTT_G[pi];

        i64 *fa = (i64 *)calloc(n_ntt, sizeof(i64));
        i64 *fb = (i64 *)calloc(n_ntt, sizeof(i64));

        for (int i = 0; i < na; i++) fa[i] = a[i] % mod;
        for (int i = 0; i < nb; i++) fb[i] = b[i] % mod;

        ntt(fa, n_ntt, 0, mod, gen);
        ntt(fb, n_ntt, 0, mod, gen);

        for (int i = 0; i < n_ntt; i++)
            fa[i] = mulmod(fa[i], fb[i], mod);

        ntt(fa, n_ntt, 1, mod, gen);
        ntt_results[pi] = fa;
        free(fb);
    }

    /* CRT via Garner's algorithm */
    for (int i = 0; i < limit; i++) {
        i64 a1 = ntt_results[0][i];  /* mod p1 */
        i64 a2 = ntt_results[1][i];  /* mod p2 */
        i64 a3 = ntt_results[2][i];  /* mod p3 */

        /* Step 1: x mod p1 = a1 */
        /* Step 2: x mod p1*p2 */
        i64 t2 = ((a2 - a1) % NTT_P[1] + NTT_P[1]) % NTT_P[1];
        t2 = mulmod(t2, p1_mod_p2_inv, NTT_P[1]);
        /* x12 = a1 + p1 * t2, but we only need x12 mod MOD */
        i64 x12_mod = (a1 % MOD + mulmod(p1_mod_MOD, t2 % MOD, MOD)) % MOD;

        /* Step 3: x mod p1*p2*p3 */
        i64 x12_mod_p3 = (a1 % NTT_P[2] + mulmod(NTT_P[0] % NTT_P[2], t2, NTT_P[2])) % NTT_P[2];
        if (x12_mod_p3 < 0) x12_mod_p3 += NTT_P[2];
        i64 t3 = ((a3 - x12_mod_p3) % NTT_P[2] + NTT_P[2]) % NTT_P[2];
        t3 = mulmod(t3, p12_mod_p3_inv, NTT_P[2]);
        /* x123 mod MOD = x12_mod + p12_mod_MOD * t3 (mod MOD) */
        res[i] = (x12_mod + mulmod(p12_mod_MOD, t3 % MOD, MOD)) % MOD;
        if (res[i] < 0) res[i] += MOD;
    }

    for (int pi = 0; pi < NTT_NPRIMES; pi++) free(ntt_results[pi]);

    return res;
}

/* ---- Series inversion via Newton iteration ---- */

static i64 *series_inverse(const i64 *f, int n_terms) {
    if (n_terms <= 0) return NULL;

    i64 f0 = f[0] % MOD;
    if (f0 < 0) f0 += MOD;

    i64 g0 = modinv(f0, MOD);
    i64 *g = (i64 *)calloc(n_terms, sizeof(i64));
    g[0] = g0;

    int m = 1;
    while (m < n_terms) {
        int m2 = m * 2;
        if (m2 > n_terms) m2 = n_terms;

        i64 *fg = conv_mod(f, m2, g, m, m2);

        for (int i = 0; i < m2; i++) {
            fg[i] = (-fg[i]) % MOD;
            if (fg[i] < 0) fg[i] += MOD;
        }
        fg[0] = (fg[0] + 2) % MOD;

        i64 *g_new = conv_mod(g, m, fg, m2, m2);

        for (int i = 0; i < m2; i++)
            g[i] = g_new[i];

        free(fg);
        free(g_new);
        m = m2;
    }

    return g;
}

/* ---- Initialize CRT coefficients ---- */

static void init_crt(void) {
    /* Garner's algorithm precomputation */
    p1_mod_p2_inv = modinv(NTT_P[0] % NTT_P[1], NTT_P[1]);
    i64 p12_mod_p3 = mulmod(NTT_P[0] % NTT_P[2], NTT_P[1] % NTT_P[2], NTT_P[2]);
    p12_mod_p3_inv = modinv(p12_mod_p3, NTT_P[2]);
    p1_mod_MOD = NTT_P[0] % MOD;
    p12_mod_MOD = mulmod(NTT_P[0] % MOD, NTT_P[1] % MOD, MOD);
}

/* ---- Solve ---- */

long long p929_native(void) {
    init_crt();

    /* Fibonacci numbers mod MOD */
    i64 *fib = (i64 *)calloc(N + 1, sizeof(i64));
    fib[1] = 1;
    if (N >= 2) fib[2] = 1;
    for (int i = 3; i <= N; i++) {
        i64 x = fib[i - 1] + fib[i - 2];
        if (x >= MOD) x -= MOD;
        fib[i] = x;
    }

    /* Sieve: s[n] = sum_{d|n} (-1)^(d-1) * Fib[d] mod MOD */
    i64 *s = (i64 *)calloc(N + 1, sizeof(i64));
    for (int d = 1; d <= N; d++) {
        i64 val = fib[d];
        if ((d & 1) == 0) val = MOD - val;
        for (int k = d; k <= N; k += d) {
            i64 x = s[k] + val;
            if (x >= MOD) x -= MOD;
            s[k] = x;
        }
    }

    /* f(x) = 1 - S(x) */
    i64 *f = (i64 *)calloc(N + 1, sizeof(i64));
    f[0] = 1;
    for (int i = 1; i <= N; i++) {
        f[i] = (MOD - s[i]) % MOD;
        if (f[i] < 0) f[i] += MOD;
    }

    /* t = 1 / f */
    i64 *t = series_inverse(f, N + 1);

    i64 result = t[N] % MOD;
    if (result < 0) result += MOD;

    free(fib); free(s); free(f); free(t);
    return result;
}
