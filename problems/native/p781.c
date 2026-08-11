// Project Euler 781: Feynman Diagrams
// Compute F(50000) mod 1e9+7 via formal power series inversion using
// three-prime NTT convolution with CRT reconstruction.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MOD 1000000007LL
#define P1  998244353LL
#define P2  1004535809LL
#define P3  469762049LL
#define G1  3LL
#define G2  3LL
#define G3  3LL

// CRT reconstruction constants (filled in p781_native).
static long long INV_P1_MOD_P2;
static long long P12_MOD_P3;
static long long INV_P12_MOD_P3;
static long long P1_MOD_P3;
static long long P1_MOD_MOD;
static long long P12_MOD_MOD;

static long long pow_mod(long long a, long long e, long long m) {
    long long r = 1, b = a % m;
    if (b < 0) b += m;
    while (e > 0) {
        if (e & 1) r = r * b % m;
        b = b * b % m;
        e >>= 1;
    }
    return r;
}

// In-place iterative NTT.
static void ntt(long long *a, int n, int invert, long long mod, long long root) {
    int lg = 0;
    while ((1 << lg) < n) lg++;
    for (int i = 1; i < n; i++) {
        int j = 0, x = i;
        for (int b = 0; b < lg; b++) { j = (j << 1) | (x & 1); x >>= 1; }
        if (i < j) { long long t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (int length = 2; length <= n; length <<= 1) {
        long long wlen = pow_mod(root, (mod - 1) / length, mod);
        if (invert) wlen = pow_mod(wlen, mod - 2, mod);
        int half = length >> 1;
        for (int i = 0; i < n; i += length) {
            long long w = 1;
            for (int j = 0; j < half; j++) {
                long long u = a[i + j];
                long long v = a[i + j + half] * w % mod;
                long long x = u + v; if (x >= mod) x -= mod;
                long long y = u - v; if (y < 0) y += mod;
                a[i + j] = x;
                a[i + j + half] = y;
                w = w * wlen % mod;
            }
        }
    }
    if (invert) {
        long long n_inv = pow_mod(n, mod - 2, mod);
        for (int i = 0; i < n; i++) a[i] = a[i] * n_inv % mod;
    }
}

// Convolution under a single NTT prime modulus. Returns malloc'd buffer.
static long long *convolution_prime(const long long *a, int an,
                                    const long long *b, int bn,
                                    long long mod, long long root, int *out_n) {
    if (an == 0 || bn == 0) { *out_n = 0; return NULL; }
    int need = an + bn - 1;
    int n = 1;
    while (n < need) n <<= 1;
    long long *fa = calloc(n, sizeof(long long));
    long long *fb = calloc(n, sizeof(long long));
    for (int i = 0; i < an; i++) fa[i] = a[i] % mod;
    for (int i = 0; i < bn; i++) fb[i] = b[i] % mod;
    ntt(fa, n, 0, mod, root);
    ntt(fb, n, 0, mod, root);
    for (int i = 0; i < n; i++) fa[i] = fa[i] * fb[i] % mod;
    ntt(fa, n, 1, mod, root);
    free(fb);
    *out_n = need;
    return fa;
}

static long long crt3_to_mod(long long r1, long long r2, long long r3) {
    long long t2 = ((r2 - r1 % P2) % P2 + P2) % P2;
    t2 = t2 * INV_P1_MOD_P2 % P2;
    long long x12_mod_p3 = (r1 % P3 + P1_MOD_P3 * t2) % P3;
    long long t3 = ((r3 - x12_mod_p3) % P3 + P3) % P3;
    t3 = t3 * INV_P12_MOD_P3 % P3;
    long long res = (r1 % MOD
                     + P1_MOD_MOD * (t2 % MOD) % MOD
                     + P12_MOD_MOD * (t3 % MOD) % MOD) % MOD;
    return res;
}

// Multiply polynomials mod MOD, writing first `limit` coefficients to res.
static void poly_mul(const long long *a, int an,
                     const long long *b, int bn,
                     int limit, long long *res) {
    if (an == 0 || bn == 0) return;
    int need = an + bn - 1;
    if (need > limit) need = limit;
    int ause = an < need ? an : need;
    int buse = bn < need ? bn : need;

    long long *a1 = malloc(ause * sizeof(long long));
    long long *b1 = malloc(buse * sizeof(long long));
    for (int i = 0; i < ause; i++) a1[i] = a[i] % P1;
    for (int i = 0; i < buse; i++) b1[i] = b[i] % P1;
    int n1; long long *c1 = convolution_prime(a1, ause, b1, buse, P1, G1, &n1);
    free(a1); free(b1);

    long long *a2 = malloc(ause * sizeof(long long));
    long long *b2 = malloc(buse * sizeof(long long));
    for (int i = 0; i < ause; i++) a2[i] = a[i] % P2;
    for (int i = 0; i < buse; i++) b2[i] = b[i] % P2;
    int n2; long long *c2 = convolution_prime(a2, ause, b2, buse, P2, G2, &n2);
    free(a2); free(b2);

    long long *a3 = malloc(ause * sizeof(long long));
    long long *b3 = malloc(buse * sizeof(long long));
    for (int i = 0; i < ause; i++) a3[i] = a[i] % P3;
    for (int i = 0; i < buse; i++) b3[i] = b[i] % P3;
    int n3; long long *c3 = convolution_prime(a3, ause, b3, buse, P3, G3, &n3);
    free(a3); free(b3);

    for (int i = 0; i < need; i++)
        res[i] = crt3_to_mod(c1[i], c2[i], c3[i]);

    free(c1); free(c2); free(c3);
}

// Formal power series inverse of a (a[0] != 0) modulo x^n.
static void poly_inv(const long long *a, int n, long long *inv) {
    inv[0] = pow_mod(a[0], MOD - 2, MOD);
    int m = 1;
    long long *t = malloc(2 * n * sizeof(long long));
    long long *u = malloc(2 * n * sizeof(long long));
    long long *newinv = malloc(2 * n * sizeof(long long));
    while (m < n) {
        int m2 = 2 * m < n ? 2 * m : n;
        poly_mul(a, m2, inv, m, m2, t);
        u[0] = (2 - t[0]) % MOD; if (u[0] < 0) u[0] += MOD;
        for (int i = 1; i < m2; i++) u[i] = (MOD - t[i]) % MOD;
        poly_mul(inv, m, u, m2, m2, newinv);
        for (int i = 0; i < m2; i++) inv[i] = newinv[i];
        m = m2;
    }
    free(t); free(u); free(newinv);
}

// Build A_m = (2m-1)!! * [x^(2m)] e^{-x}/(1-x) and B_m with /(1-x)^2.
static void build_series(int M, long long *A, long long *B) {
    int N = 2 * M;
    long long *fact = malloc((N + 1) * sizeof(long long));
    long long *invfact = malloc((N + 1) * sizeof(long long));
    long long *S = malloc((N + 1) * sizeof(long long));

    fact[0] = 1;
    for (int i = 1; i <= N; i++) fact[i] = fact[i - 1] * i % MOD;
    invfact[N] = pow_mod(fact[N], MOD - 2, MOD);
    for (int i = N; i >= 1; i--) invfact[i - 1] = invfact[i] * i % MOD;

    long long acc = 0;
    for (int r = 0; r <= N; r++) {
        long long term = invfact[r];
        if (r & 1) term = MOD - term;
        acc += term; if (acc >= MOD) acc -= MOD;
        S[r] = acc;
    }

    long long df = 1;
    for (int m = 0; m <= M; m++) {
        if (m > 0) df = df * (2 * m - 1) % MOD;
        long long s2m = S[2 * m];
        long long s2m_1 = (2 * m - 1 >= 0) ? S[2 * m - 1] : 0;
        long long a = s2m;
        long long b = ((2 * m + 1) * s2m + s2m_1) % MOD;
        A[m] = df * a % MOD;
        B[m] = df * b % MOD;
    }
    free(fact); free(invfact); free(S);
}

long long p781_native(void) {
    // CRT constants.
    INV_P1_MOD_P2 = pow_mod(P1, P2 - 2, P2);
    P12_MOD_P3 = (P1 * P2) % P3;
    INV_P12_MOD_P3 = pow_mod(P12_MOD_P3, P3 - 2, P3);
    P1_MOD_P3 = P1 % P3;
    P1_MOD_MOD = P1 % MOD;
    P12_MOD_MOD = (P1 % MOD) * (P2 % MOD) % MOD;

    int n = 50000;
    int M = n / 2;

    long long *A = malloc((M + 1) * sizeof(long long));
    long long *B = malloc((M + 1) * sizeof(long long));
    build_series(M, A, B);

    long long *invA = malloc((M + 1) * sizeof(long long));
    poly_inv(A, M + 1, invA);

    long long *G = malloc((M + 1) * sizeof(long long));
    poly_mul(B, M + 1, invA, M + 1, M + 1, G);

    long long ans = G[M] % MOD;

    free(A); free(B); free(invA); free(G);
    return ans;
}
