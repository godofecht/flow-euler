
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

#define MOD 1000000123LL

static i64 *poly_mul_mod(const i64 *A, int n, const i64 *B, int m) {
    /* naive for small */
    if ((i64)n * m <= 40000 || (n < 40) || (m < 40)) {
        i64 *res = calloc((size_t)(n + m - 1), sizeof(i64));
        for (int i = 0; i < n; i++) if (A[i]) {
            for (int j = 0; j < m; j++) {
                res[i + j] = (res[i + j] + (i64)((__int128)A[i] * B[j] % MOD)) % MOD;
            }
        }
        return res;
    }
    int on;
    return convolve_mod(A, n, B, m, MOD, &on);
}

static i64 *poly_inv_mod(const i64 *f, int n) {
    i64 *g = calloc((size_t)n, sizeof(i64));
    g[0] = mod_inv(f[0], MOD);
    int m = 1;
    while (m < n) {
        int m2 = m * 2; if (m2 > n) m2 = n;
        i64 *fg = poly_mul_mod(f, m2, g, m);
        for (int i = 0; i < m2; i++) {
            i64 v = (i < m + m - 1) ? fg[i] : 0; /* length may be m2+m-1 */
            (void)v;
        }
        /* trim fg to m2 */
        i64 *h = calloc((size_t)m2, sizeof(i64));
        int fglen = m2 + m - 1;
        for (int i = 0; i < m2; i++) {
            i64 v = (i < fglen) ? fg[i] : 0;
            h[i] = (MOD - v) % MOD;
        }
        h[0] = (h[0] + 2) % MOD;
        free(fg);
        i64 *ng = poly_mul_mod(g, m, h, m2);
        free(g); free(h);
        g = calloc((size_t)m2, sizeof(i64));
        int nglen = m + m2 - 1;
        for (int i = 0; i < m2; i++) g[i] = (i < nglen) ? ng[i] : 0;
        free(ng);
        m = m2;
    }
    return g;
}

static void pre_fac_pow(int n, int r, i64 **fac_pow, i64 **inv_fac_pow) {
    *fac_pow = malloc((size_t)(n + 1) * sizeof(i64));
    *inv_fac_pow = malloc((size_t)(n + 1) * sizeof(i64));
    (*fac_pow)[0] = 1;
    for (int i = 1; i <= n; i++) {
        i64 ip = mod_pow(i, r, MOD);
        (*fac_pow)[i] = (i64)((__int128)(*fac_pow)[i - 1] * ip % MOD);
    }
    (*inv_fac_pow)[n] = mod_inv((*fac_pow)[n], MOD);
    for (int i = n - 1; i >= 0; i--) {
        i64 ip = mod_pow(i + 1, r, MOD);
        (*inv_fac_pow)[i] = (i64)((__int128)(*inv_fac_pow)[i + 1] * ip % MOD);
    }
}

static i64 P_via_dp(int k, int n, i64 *fac_pow, i64 *inv_fac_pow) {
    int q = n / k, rem = n % k;
    int m = q + (rem ? 1 : 0);
    int *blocks = malloc((size_t)m * sizeof(int));
    for (int i = 0; i < q; i++) blocks[i] = k;
    if (rem) blocks[q] = rem;
    i64 *f = calloc((size_t)(m + 1), sizeof(i64));
    f[0] = fac_pow[n];
    for (int i = 1; i <= m; i++) {
        int seg = 0;
        i64 acc = 0;
        for (int j = i - 1; j >= 0; j--) {
            seg += blocks[j];
            i64 term = (i64)((__int128)f[j] * inv_fac_pow[seg] % MOD);
            if ((i - j) & 1) acc += term; else acc -= term;
        }
        acc %= MOD; if (acc < 0) acc += MOD;
        f[i] = acc;
    }
    i64 ans = f[m];
    free(blocks); free(f);
    return ans;
}

static i64 P_via_inverse(int k, int n, i64 *fac_pow, i64 *inv_fac_pow) {
    int q = n / k, rem = n % k;
    i64 *p = calloc((size_t)(q + 1), sizeof(i64));
    p[0] = 1;
    for (int d = 1; d <= q; d++) {
        i64 coeff = inv_fac_pow[d * k];
        if (d & 1) coeff = (MOD - coeff) % MOD;
        p[d] = coeff;
    }
    i64 *s = poly_inv_mod(p, q + 1);
    i64 g0 = fac_pow[n];
    i64 *g = malloc((size_t)(q + 1) * sizeof(i64));
    for (int i = 0; i <= q; i++) g[i] = (i64)((__int128)g0 * s[i] % MOD);
    i64 ans;
    if (rem == 0) ans = g[q];
    else {
        ans = 0;
        for (int t = 0; t <= q; t++) {
            int seg_len = (q - t) * k + rem;
            i64 term = (i64)((__int128)g[t] * inv_fac_pow[seg_len] % MOD);
            if ((q - t) & 1) ans -= term; else ans += term;
        }
        ans %= MOD; if (ans < 0) ans += MOD;
    }
    free(p); free(s); free(g);
    return ans;
}

static i64 Qn(int n) {
    int r = n;
    i64 *fac_pow, *inv_fac_pow;
    pre_fac_pow(n, r, &fac_pow, &inv_fac_pow);
    int k_inv_max = n / 200;
    i64 ans = 0;
    for (int k = 1; k <= n; k++) {
        int q = n / k;
        i64 v;
        if (k <= k_inv_max && q > 200) v = P_via_inverse(k, n, fac_pow, inv_fac_pow);
        else v = P_via_dp(k, n, fac_pow, inv_fac_pow);
        ans += v; if (ans >= MOD) ans -= MOD;
    }
    free(fac_pow); free(inv_fac_pow);
    return ans;
}

long long pe559_answer(void) {
    return Qn(50000);
}
