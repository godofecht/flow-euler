/* Project Euler 840 - Sum of Products.
 * Computes S(5*10^4) mod 999676999 using CDQ divide-and-conquer with NTT.
 *
 * G(0)=1, G(k) = (1/k) * sum_{j=1..k} b[j]*G(k-j)
 * b[k] = sum_{m|k} m * D(m)^{k/m}  where D is the arithmetic derivative.
 * Answer = sum_{k=1..N} G(k) mod MOD.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;
typedef __int128 i128;

#define MOD 999676999LL
#define N 50000
#define THRESH 256

static const i64 P1 = 998244353, P2 = 1004535809, P3 = 469762049;

/* ---------- NTT ---------- */

static i64 mod_pow(i64 a, i64 e, i64 mod) {
    i64 r = 1 % mod; a %= mod;
    while (e > 0) {
        if (e & 1) r = (i64)((i128)r * a % mod);
        a = (i64)((i128)a * a % mod);
        e >>= 1;
    }
    return r;
}

static i64 mod_inv(i64 a, i64 mod) { return mod_pow(a, mod - 2, mod); }

static void ntt(i64 *a, int n, int invert, i64 mod, i64 root) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { i64 t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        i64 wlen = mod_pow(root, (mod - 1) / len, mod);
        if (invert) wlen = mod_inv(wlen, mod);
        for (int i = 0; i < n; i += len) {
            i64 w = 1;
            for (int j = 0; j < len / 2; j++) {
                i64 u = a[i + j];
                i64 v = (i64)((i128)a[i + j + len / 2] * w % mod);
                i64 x = u + v; if (x >= mod) x -= mod;
                i64 y = u - v; if (y < 0) y += mod;
                a[i + j] = x; a[i + j + len / 2] = y;
                w = (i64)((i128)w * wlen % mod);
            }
        }
    }
    if (invert) {
        i64 ninv = mod_inv(n, mod);
        for (int i = 0; i < n; i++) a[i] = (i64)((i128)a[i] * ninv % mod);
    }
}

/* ---------- CRT ---------- */

static i64 inv_p1_p2, inv_p12_p3, p1_mod, p12_mod;

static void crt_init(void) {
    inv_p1_p2 = mod_inv(P1 % P2, P2);
    i128 p12 = (i128)P1 * P2;
    inv_p12_p3 = mod_inv((i64)(p12 % P3), P3);
    p1_mod = P1 % MOD;
    p12_mod = (i64)(p12 % MOD);
}

static i64 crt3(i64 a1, i64 a2, i64 a3) {
    i64 t1 = (a2 - a1) % P2; if (t1 < 0) t1 += P2;
    t1 = (i64)((i128)t1 * inv_p1_p2 % P2);
    i64 x12_p3 = (i64)(((i128)a1 + (i128)P1 * t1) % P3);
    i64 t2 = (a3 - x12_p3) % P3; if (t2 < 0) t2 += P3;
    t2 = (i64)((i128)t2 * inv_p12_p3 % P3);
    i128 r = (i128)(a1 % MOD) + (i128)p1_mod * (t1 % MOD) + (i128)p12_mod * (t2 % MOD);
    return (i64)(r % MOD);
}

/* ---------- NTT cache for b ---------- */

typedef struct { int L; i64 *fb[3]; } NTTCacheEntry;
static NTTCacheEntry ntt_cache[20];
static int ntt_cache_count = 0;

static i64 *g_b;  /* b array */

static NTTCacheEntry *get_ntt_cache(int L) {
    for (int i = 0; i < ntt_cache_count; i++)
        if (ntt_cache[i].L == L) return &ntt_cache[i];
    NTTCacheEntry *c = &ntt_cache[ntt_cache_count++];
    c->L = L;
    int nfft = 2 * L;
    i64 primes[3] = {P1, P2, P3};
    for (int pi = 0; pi < 3; pi++) {
        c->fb[pi] = (i64 *)calloc((size_t)nfft, sizeof(i64));
        for (int i = 0; i < L; i++) c->fb[pi][i] = g_b[i] % primes[pi];
        ntt(c->fb[pi], nfft, 0, primes[pi], 3);
    }
    return c;
}

/* Convolve a_seg[0..half-1] with b[0..L-1], return first L coefficients mod MOD. */
static void convolve_first_L(i64 *a_seg, int half, int L, i64 *out) {
    int nfft = 2 * L;
    i64 primes[3] = {P1, P2, P3};
    NTTCacheEntry *c = get_ntt_cache(L);

    i64 *fa = (i64 *)calloc((size_t)nfft, sizeof(i64));
    i64 *residues[3];
    for (int pi = 0; pi < 3; pi++)
        residues[pi] = (i64 *)malloc((size_t)L * sizeof(i64));

    for (int pi = 0; pi < 3; pi++) {
        i64 p = primes[pi];
        memset(fa, 0, (size_t)nfft * sizeof(i64));
        for (int i = 0; i < half; i++) fa[i] = a_seg[i] % p;
        ntt(fa, nfft, 0, p, 3);
        for (int i = 0; i < nfft; i++)
            fa[i] = (i64)((i128)fa[i] * c->fb[pi][i] % p);
        ntt(fa, nfft, 1, p, 3);
        for (int i = 0; i < L; i++) residues[pi][i] = fa[i];
    }

    for (int i = 0; i < L; i++)
        out[i] = crt3(residues[0][i], residues[1][i], residues[2][i]);

    free(fa);
    for (int pi = 0; pi < 3; pi++) free(residues[pi]);
}

/* ---------- Number theory ---------- */

static int next_pow2(int x) { int p = 1; while (p < x) p <<= 1; return p; }

static void sieve_spf(int *spf, int limit) {
    for (int i = 0; i <= limit; i++) spf[i] = i;
    spf[0] = 0; if (limit >= 1) spf[1] = 1;
    for (int i = 2; (i64)i * i <= limit; i++)
        if (spf[i] == i)
            for (int j = i * i; j <= limit; j += i)
                if (spf[j] == j) spf[j] = i;
}

/* ---------- Global solver state ---------- */

static i64 *g_a, *g_f, *g_inv;
static int g_n, g_size;

static void solve_block(int l, int r) {
    int rr = r < g_n ? r : g_n;
    for (int i = l; i <= rr; i++) {
        i64 ai = (i == 0) ? 1 : (i64)((i128)g_f[i] * g_inv[i] % MOD);
        g_a[i] = ai;
        int max_k = rr - i;
        if (max_k <= 0 || ai == 0) continue;
        for (int k = 1; k <= max_k; k++)
            g_f[i + k] = (g_f[i + k] + (i64)((i128)ai * g_b[k] % MOD)) % MOD;
    }
}

static void cdq(int l, int L) {
    if (l > g_n) return;
    int r = l + L - 1;
    if (L <= THRESH) { solve_block(l, r); return; }
    int half = L >> 1;
    int mid = l + half - 1;

    cdq(l, half);

    if (mid < g_n) {
        i64 *conv = (i64 *)malloc((size_t)L * sizeof(i64));
        convolve_first_L(&g_a[l], half, L, conv);
        int rr = r < g_n ? r : g_n;
        for (int t = mid + 1; t <= rr; t++)
            g_f[t] = (g_f[t] + conv[t - l]) % MOD;
        free(conv);
    }

    cdq(mid + 1, half);
}

long long p840_native(void) {
    crt_init();
    g_n = N;
    g_size = next_pow2(N + 1);

    /* SPF and arithmetic derivative D */
    int *spf = (int *)malloc((size_t)(N + 1) * sizeof(int));
    sieve_spf(spf, N);
    i64 *D = (i64 *)calloc((size_t)(N + 1), sizeof(i64));
    for (int n = 2; n <= N; n++) {
        int p = spf[n], m = n / p;
        D[n] = m + (i64)p * D[m];
    }
    D[1] = 1;

    /* Build b[k] = sum_{m|k} m * D(m)^{k/m} mod MOD */
    g_b = (i64 *)calloc((size_t)g_size, sizeof(i64));
    for (int m = 1; m <= N; m++) {
        i64 wm = D[m] % MOD;
        i64 mm = m % MOD;
        i64 pwr = wm;
        for (int j = m; j <= N; j += m) {
            g_b[j] = (g_b[j] + (i64)((i128)mm * pwr % MOD)) % MOD;
            pwr = (i64)((i128)pwr * wm % MOD);
        }
    }

    /* Modular inverses */
    g_inv = (i64 *)calloc((size_t)(N + 1), sizeof(i64));
    g_inv[1] = 1;
    for (int i = 2; i <= N; i++)
        g_inv[i] = (MOD - (MOD / i) * g_inv[MOD % i] % MOD) % MOD;

    /* CDQ */
    g_a = (i64 *)calloc((size_t)g_size, sizeof(i64));
    g_f = (i64 *)calloc((size_t)g_size, sizeof(i64));
    cdq(0, g_size);

    /* Answer = sum(G[1..N]) mod MOD */
    i64 ans = 0;
    for (int k = 1; k <= N; k++)
        ans = (ans + g_a[k]) % MOD;

    /* Cleanup */
    for (int i = 0; i < ntt_cache_count; i++)
        for (int j = 0; j < 3; j++) free(ntt_cache[i].fb[j]);
    ntt_cache_count = 0;
    free(spf); free(D); free(g_b); free(g_inv); free(g_a); free(g_f);

    return (long long)ans;
}
