
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef int64_t i64;

static i64 mod_pow(i64 a, i64 e, i64 mod) {
    i64 r = 1 % mod;
    a %= mod;
    while (e > 0) {
        if (e & 1) r = (i64)((__int128)r * a % mod);
        a = (i64)((__int128)a * a % mod);
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
                i64 v = (i64)((__int128)a[i + j + len / 2] * w % mod);
                i64 x = u + v; if (x >= mod) x -= mod;
                i64 y = u - v; if (y < 0) y += mod;
                a[i + j] = x;
                a[i + j + len / 2] = y;
                w = (i64)((__int128)w * wlen % mod);
            }
        }
    }
    if (invert) {
        i64 ninv = mod_inv(n, mod);
        for (int i = 0; i < n; i++) a[i] = (i64)((__int128)a[i] * ninv % mod);
    }
}

static i64 *convolve_one(const i64 *a, int n1, const i64 *b, int n2, i64 mod, i64 root, int *out_n) {
    int need = n1 + n2 - 1;
    int n = 1; while (n < need) n <<= 1;
    i64 *fa = calloc((size_t)n, sizeof(i64));
    i64 *fb = calloc((size_t)n, sizeof(i64));
    for (int i = 0; i < n1; i++) fa[i] = a[i] % mod;
    for (int i = 0; i < n2; i++) fb[i] = b[i] % mod;
    ntt(fa, n, 0, mod, root);
    ntt(fb, n, 0, mod, root);
    for (int i = 0; i < n; i++) fa[i] = (i64)((__int128)fa[i] * fb[i] % mod);
    ntt(fa, n, 1, mod, root);
    *out_n = need;
    i64 *res = malloc((size_t)need * sizeof(i64));
    memcpy(res, fa, (size_t)need * sizeof(i64));
    free(fa); free(fb);
    return res;
}

/* CRT combine three NTT moduli into result mod MOD */
static const i64 P1 = 998244353, G1 = 3;
static const i64 P2 = 1004535809, G2 = 3;
static const i64 P3 = 469762049, G3 = 3;

static i64 crt3(i64 a1, i64 a2, i64 a3, i64 MOD) {
    static i64 inv_p1_mod_p2 = -1, inv_p12_mod_p3 = -1;
    static __int128 P12;
    if (inv_p1_mod_p2 < 0) {
        inv_p1_mod_p2 = mod_inv(P1 % P2, P2);
        P12 = (__int128)P1 * P2;
        inv_p12_mod_p3 = mod_inv((i64)(P12 % P3), P3);
    }
    i64 t1 = (i64)((__int128)(a2 - a1) % P2);
    if (t1 < 0) t1 += P2;
    t1 = (i64)((__int128)t1 * inv_p1_mod_p2 % P2);
    __int128 x2 = a1 + (__int128)P1 * t1;
    i64 t2 = (i64)((a3 - (i64)(x2 % P3)) % P3);
    if (t2 < 0) t2 += P3;
    t2 = (i64)((__int128)t2 * inv_p12_mod_p3 % P3);
    __int128 x = x2 + P12 * t2;
    i64 r = (i64)(x % MOD);
    if (r < 0) r += MOD;
    return r;
}

static i64 *convolve_mod(const i64 *a, int n1, const i64 *b, int n2, i64 MOD, int *out_n) {
    if (n1 == 0 || n2 == 0) { *out_n = 0; return NULL; }
    int o1,o2,o3;
    i64 *c1 = convolve_one(a,n1,b,n2,P1,G1,&o1);
    i64 *c2 = convolve_one(a,n1,b,n2,P2,G2,&o2);
    i64 *c3 = convolve_one(a,n1,b,n2,P3,G3,&o3);
    int need = n1 + n2 - 1;
    i64 *res = malloc((size_t)need * sizeof(i64));
    for (int i = 0; i < need; i++) res[i] = crt3(c1[i], c2[i], c3[i], MOD);
    free(c1); free(c2); free(c3);
    *out_n = need;
    return res;
}

#define MOD 1000000007LL

static i64 *poly_inv(const i64 *f, int deg) {
    i64 *g = calloc((size_t)deg, sizeof(i64));
    g[0] = mod_inv(f[0], MOD);
    int m = 1;
    while (m < deg) {
        int m2 = m * 2; if (m2 > deg) m2 = deg;
        int on;
        i64 *fg = convolve_mod(f, m2, g, m, MOD, &on);
        if (on > m2) on = m2;
        i64 *h = calloc((size_t)m2, sizeof(i64));
        for (int i = 0; i < on && i < m2; i++) h[i] = (MOD - fg[i]) % MOD;
        h[0] = (h[0] + 2) % MOD;
        free(fg);
        i64 *ng = convolve_mod(g, m, h, m2, MOD, &on);
        free(g); free(h);
        g = calloc((size_t)m2, sizeof(i64));
        for (int i = 0; i < m2 && i < on; i++) g[i] = ng[i];
        free(ng);
        m = m2;
    }
    return g;
}

static i64 *poly_log(const i64 *f, int deg, const i64 *inv_int) {
    /* f[0]==1; res[0]=0; res[i] = (f'/f)[i-1] / i */
    i64 *df = malloc((size_t)(deg - 1) * sizeof(i64));
    int dflen = 0;
    for (int i = 1; i < deg; i++) {
        df[i - 1] = (i64)((__int128)i * f[i] % MOD);
        dflen++;
    }
    i64 *invf = poly_inv(f, deg);
    int on;
    i64 *prod = convolve_mod(df, dflen, invf, deg, MOD, &on);
    i64 *res = calloc((size_t)deg, sizeof(i64));
    for (int i = 1; i < deg; i++) {
        if (i - 1 < on) res[i] = (i64)((__int128)prod[i - 1] * inv_int[i] % MOD);
    }
    free(df); free(invf); free(prod);
    return res;
}

static i64 *poly_pow(const i64 *a, int deg, int exp) {
    i64 *res = calloc((size_t)deg, sizeof(i64));
    res[0] = 1;
    i64 *base = malloc((size_t)deg * sizeof(i64));
    memcpy(base, a, (size_t)deg * sizeof(i64));
    int e = exp;
    while (e > 0) {
        if (e & 1) {
            int on;
            i64 *tmp = convolve_mod(res, deg, base, deg, MOD, &on);
            memset(res, 0, (size_t)deg * sizeof(i64));
            for (int i = 0; i < deg && i < on; i++) res[i] = tmp[i];
            free(tmp);
        }
        e >>= 1;
        if (e) {
            int on;
            i64 *tmp = convolve_mod(base, deg, base, deg, MOD, &on);
            memset(base, 0, (size_t)deg * sizeof(i64));
            for (int i = 0; i < deg && i < on; i++) base[i] = tmp[i];
            free(tmp);
        }
    }
    free(base);
    return res;
}

static void precompute(int n, i64 **fact, i64 **inv_fact, i64 **inv_int) {
    *fact = malloc((size_t)(n + 1) * sizeof(i64));
    *inv_fact = malloc((size_t)(n + 1) * sizeof(i64));
    *inv_int = malloc((size_t)(n + 1) * sizeof(i64));
    (*fact)[0] = 1;
    for (int i = 1; i <= n; i++) (*fact)[i] = (i64)((__int128)(*fact)[i - 1] * i % MOD);
    (*inv_fact)[n] = mod_inv((*fact)[n], MOD);
    for (int i = n; i > 0; i--) (*inv_fact)[i - 1] = (i64)((__int128)(*inv_fact)[i] * i % MOD);
    (*inv_int)[0] = 0;
    if (n >= 1) (*inv_int)[1] = 1;
    for (int i = 2; i <= n; i++) (*inv_int)[i] = MOD - (i64)((__int128)(MOD / i) * (*inv_int)[MOD % i] % MOD);
}

static i64 *compute_A(int max_n, i64 *inv_fact, i64 *inv_int) {
    i64 modm1 = MOD - 1;
    i64 *exp2 = malloc((size_t)(max_n + 1) * sizeof(i64));
    exp2[0] = 1;
    for (int i = 1; i <= max_n; i++) exp2[i] = (exp2[i - 1] * 2) % modm1;
    i64 *A0 = malloc((size_t)(max_n + 1) * sizeof(i64));
    for (int i = 0; i <= max_n; i++) A0[i] = mod_pow(2, (exp2[i] - 1) % modm1, MOD);
    i64 *p = malloc((size_t)(max_n + 1) * sizeof(i64));
    i64 *q = malloc((size_t)(max_n + 1) * sizeof(i64));
    for (int i = 0; i <= max_n; i++) {
        p[i] = (i64)((__int128)A0[i] * inv_fact[i] % MOD);
        q[i] = (i & 1) ? (MOD - inv_fact[i]) % MOD : inv_fact[i];
    }
    int on;
    i64 *H = convolve_mod(p, max_n + 1, q, max_n + 1, MOD, &on);
    i64 *H2 = calloc((size_t)(max_n + 1), sizeof(i64));
    for (int i = 0; i <= max_n && i < on; i++) H2[i] = H[i];
    H2[0] = 1;
    free(H);
    i64 *A = poly_log(H2, max_n + 1, inv_int);
    free(exp2); free(A0); free(p); free(q); free(H2);
    return A;
}

static i64 C_nk(int n, int k, i64 *A_series, i64 *fact, i64 *inv_fact) {
    int deg = n + 1;
    i64 *Ak = poly_pow(A_series, deg, k);
    int on;
    i64 *prod = convolve_mod(Ak, deg, inv_fact, deg, MOD, &on);
    i64 ans = (i64)((__int128)fact[n] * (n < on ? prod[n] : 0) % MOD * inv_fact[k] % MOD);
    free(Ak); free(prod);
    return ans;
}

long long pe553_answer(void) {
    int N = 10000, K = 10;
    i64 *fact, *inv_fact, *inv_int;
    precompute(N, &fact, &inv_fact, &inv_int);
    i64 *A = compute_A(N, inv_fact, inv_int);
    i64 ans = C_nk(N, K, A, fact, inv_fact);
    free(fact); free(inv_fact); free(inv_int); free(A);
    return ans;
}
