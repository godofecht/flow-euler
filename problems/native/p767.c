#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;
typedef __int128 i128;

enum { MOD = 1000000007LL };

/* Three NTT-friendly primes for exact integer convolution reconstruction. */
static const i64 PRIMES[3] = { 998244353LL, 1004535809LL, 469762049LL };
static const i64 PRIMES_G[3] = { 3, 3, 3 };

static i64 ceil_pow2(i64 x) {
    i64 n = 1;
    while (n < x) n <<= 1;
    return n;
}

static i64 mod_pow(i64 a, i64 e, i64 mod) {
    i64 r = 1; a %= mod; if (a < 0) a += mod;
    while (e > 0) {
        if (e & 1) r = (i64)((i128)r * a % mod);
        a = (i64)((i128)a * a % mod);
        e >>= 1;
    }
    return r;
}

/* ---- NTT context ---- */
typedef struct {
    i64 n;
    int *rev;
    i64 *roots[3];
    i64 *roots_inv[3];
    /* CRT constants */
    i64 inv_p1_mod_p2;
    i64 p1_mod_p3, p12_mod_p3, inv_p12_mod_p3;
    i64 p1_mod_mod, p12_mod_mod;
} NTTCtx;

static void ntt(i64 *a, i64 n, i64 mod, const i64 *roots, const i64 *roots_inv,
                const int *rev, int invert) {
    for (i64 i = 0; i < n; i++) {
        int j = rev[i];
        if (i < (i64)j) { i64 t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    const i64 *r = invert ? roots_inv : roots;
    for (i64 length = 2; length <= n; length <<= 1) {
        i64 half = length >> 1;
        i64 step = n / length;
        for (i64 i0 = 0; i0 < n; i0 += length) {
            i64 idx = 0;
            for (i64 j = i0; j < i0 + half; j++) {
                i64 u = a[j];
                i64 v = (i64)((i128)a[j + half] * r[idx] % mod);
                i64 x = u + v; if (x >= mod) x -= mod;
                i64 y = u - v; if (y < 0) y += mod;
                a[j] = x;
                a[j + half] = y;
                idx += step;
            }
        }
    }
    if (invert) {
        i64 n_inv = mod_pow(n, mod - 2, mod);
        for (i64 i = 0; i < n; i++)
            a[i] = (i64)((i128)a[i] * n_inv % mod);
    }
}

static void ctx_init(NTTCtx *c, i64 n) {
    c->n = n;
    c->rev = malloc(sizeof(int) * n);
    {
        int j = 0;
        for (i64 i = 1; i < n; i++) {
            int bit = (int)(n >> 1);
            while (j & bit) { j ^= bit; bit >>= 1; }
            j ^= bit;
            c->rev[i] = j;
        }
        c->rev[0] = 0;
    }
    for (int t = 0; t < 3; t++) {
        i64 p = PRIMES[t], g = PRIMES_G[t];
        i64 w = mod_pow(g, (p - 1) / n, p);
        i64 w_inv = mod_pow(w, p - 2, p);
        c->roots[t] = malloc(sizeof(i64) * n);
        c->roots_inv[t] = malloc(sizeof(i64) * n);
        c->roots[t][0] = 1;
        c->roots_inv[t][0] = 1;
        for (i64 i = 1; i < n; i++) {
            c->roots[t][i] = (i64)((i128)c->roots[t][i-1] * w % p);
            c->roots_inv[t][i] = (i64)((i128)c->roots_inv[t][i-1] * w_inv % p);
        }
    }
    i64 p1 = PRIMES[0], p2 = PRIMES[1], p3 = PRIMES[2];
    c->inv_p1_mod_p2 = mod_pow(p1, p2 - 2, p2);
    c->p1_mod_p3 = p1 % p3;
    c->p12_mod_p3 = (i64)((i128)p1 * p2 % p3);
    c->inv_p12_mod_p3 = mod_pow(c->p12_mod_p3, p3 - 2, p3);
    c->p1_mod_mod = p1 % MOD;
    c->p12_mod_mod = (i64)((i128)p1 * p2 % MOD);
}

static void ctx_free(NTTCtx *c) {
    free(c->rev);
    for (int t = 0; t < 3; t++) { free(c->roots[t]); free(c->roots_inv[t]); }
}

/* CRT-reduce three residue arrays (length need) to values mod MOD. */
static void crt_reduce(NTTCtx *c, i64 *out, i64 *convs[3], i64 need) {
    i64 p1 = PRIMES[0], p2 = PRIMES[1], p3 = PRIMES[2];
    for (i64 i = 0; i < need; i++) {
        i64 r1 = convs[0][i];
        i64 r2 = convs[1][i];
        i64 r3 = convs[2][i];
        i64 t1 = ((r2 - r1) % p2 + p2) % p2;
        t1 = (i64)((i128)t1 * c->inv_p1_mod_p2 % p2);
        i64 x2_mod_p3 = (r1 + (i64)((i128)c->p1_mod_p3 * t1 % p3)) % p3;
        i64 t2 = ((r3 - x2_mod_p3) % p3 + p3) % p3;
        t2 = (i64)((i128)t2 * c->inv_p12_mod_p3 % p3);
        i64 v = (r1 + (i64)((i128)c->p1_mod_mod * (t1 % MOD) % MOD)
                    + (i64)((i128)c->p12_mod_mod * (t2 % MOD) % MOD)) % MOD;
        out[i] = v;
    }
    (void)p1; (void)p2; (void)p3;
}

/* Convolution of a and b (lengths la, lb) mod MOD via 3-prime NTT + CRT. */
static i64 *convolution_mod(NTTCtx *c, const i64 *a, i64 la,
                            const i64 *b, i64 lb, i64 *need_out) {
    i64 need = la + lb - 1;
    i64 n = c->n;
    i64 *convs[3];
    for (int t = 0; t < 3; t++) {
        i64 p = PRIMES[t];
        i64 *fa = calloc(n, sizeof(i64));
        i64 *fb = calloc(n, sizeof(i64));
        for (i64 i = 0; i < la; i++) fa[i] = ((a[i] % p) + p) % p;
        for (i64 i = 0; i < lb; i++) fb[i] = ((b[i] % p) + p) % p;
        ntt(fa, n, p, c->roots[t], c->roots_inv[t], c->rev, 0);
        ntt(fb, n, p, c->roots[t], c->roots_inv[t], c->rev, 0);
        for (i64 i = 0; i < n; i++)
            fa[i] = (i64)((i128)fa[i] * fb[i] % p);
        ntt(fa, n, p, c->roots[t], c->roots_inv[t], c->rev, 1);
        convs[t] = malloc(sizeof(i64) * need);
        memcpy(convs[t], fa, sizeof(i64) * need);
        free(fa); free(fb);
    }
    i64 *out = malloc(sizeof(i64) * need);
    crt_reduce(c, out, convs, need);
    for (int t = 0; t < 3; t++) free(convs[t]);
    *need_out = need;
    return out;
}

/* Square convolution: a*a. */
static i64 *convolution_square_mod(NTTCtx *c, const i64 *a, i64 la, i64 *need_out) {
    i64 need = 2 * la - 1;
    i64 n = c->n;
    i64 *convs[3];
    for (int t = 0; t < 3; t++) {
        i64 p = PRIMES[t];
        i64 *fa = calloc(n, sizeof(i64));
        for (i64 i = 0; i < la; i++) fa[i] = ((a[i] % p) + p) % p;
        ntt(fa, n, p, c->roots[t], c->roots_inv[t], c->rev, 0);
        for (i64 i = 0; i < n; i++)
            fa[i] = (i64)((i128)fa[i] * fa[i] % p);
        ntt(fa, n, p, c->roots[t], c->roots_inv[t], c->rev, 1);
        convs[t] = malloc(sizeof(i64) * need);
        memcpy(convs[t], fa, sizeof(i64) * need);
        free(fa);
    }
    i64 *out = malloc(sizeof(i64) * need);
    crt_reduce(c, out, convs, need);
    for (int t = 0; t < 3; t++) free(convs[t]);
    *need_out = need;
    return out;
}

static i64 pow16_mod(i64 x) {
    i64 x2  = (i64)((i128)x * x % MOD);
    i64 x4  = (i64)((i128)x2 * x2 % MOD);
    i64 x8  = (i64)((i128)x4 * x4 % MOD);
    return (i64)((i128)x8 * x8 % MOD);
}

long long p767_native(void) {
    i64 k = 100000;
    i64 n = (i64)10000000000000000LL; /* 10^16 */
    /* n is divisible by k: m = 10^11 */
    i64 m = n / k;
    i64 A = mod_pow(2, m % (MOD - 1), MOD); /* 2^m mod MOD (Fermat) */

    i64 *fact = malloc(sizeof(i64) * (k + 1));
    i64 *inv_fact = malloc(sizeof(i64) * (k + 1));
    fact[0] = 1;
    for (i64 i = 1; i <= k; i++)
        fact[i] = (i64)((i128)fact[i-1] * i % MOD);
    inv_fact[k] = mod_pow(fact[k], MOD - 2, MOD);
    for (i64 i = k; i >= 1; i--)
        inv_fact[i-1] = (i64)((i128)inv_fact[i] * i % MOD);

    i64 *fact16 = malloc(sizeof(i64) * (k + 1));
    i64 *inv_fact16 = malloc(sizeof(i64) * (k + 1));
    for (i64 i = 0; i <= k; i++) {
        fact16[i] = pow16_mod(fact[i]);
        inv_fact16[i] = pow16_mod(inv_fact[i]);
    }

    i64 need = 2 * (k + 1) - 1;
    i64 ntt_len = ceil_pow2(need);
    NTTCtx ctx;
    ctx_init(&ctx, ntt_len);

    /* f[r] = sum_x C(r,x)^16 = r!^16 * [y^r] (sum_n y^n / n!^16)^2 */
    i64 conv_aa_len;
    i64 *conv_aa = convolution_square_mod(&ctx, inv_fact16, k + 1, &conv_aa_len);

    i64 *f = malloc(sizeof(i64) * (k + 1));
    for (i64 r = 0; r <= k; r++)
        f[r] = (i64)((i128)fact16[r] * conv_aa[r] % MOD);
    free(conv_aa);

    /* Binomial transform with parameter -2:
       S[L] = sum_{r<=L} C(L,r) (-2)^{L-r} f[r]
       S[L]/L! = sum_{r<=L} (f[r]/r!) * ((-2)^{L-r}/(L-r)!) */
    i64 *fprime = malloc(sizeof(i64) * (k + 1));
    for (i64 r = 0; r <= k; r++)
        fprime[r] = (i64)((i128)f[r] * inv_fact[r] % MOD);

    i64 neg2 = MOD - 2;
    i64 *b = malloc(sizeof(i64) * (k + 1));
    i64 p = 1;
    for (i64 j = 0; j <= k; j++) {
        b[j] = (i64)((i128)p * inv_fact[j] % MOD);
        p = (i64)((i128)p * neg2 % MOD);
    }

    i64 conv_fb_len;
    i64 *conv_fb = convolution_mod(&ctx, fprime, k + 1, b, k + 1, &conv_fb_len);

    i64 *S = malloc(sizeof(i64) * (k + 1));
    for (i64 L = 0; L <= k; L++)
        S[L] = (i64)((i128)conv_fb[L] * fact[L] % MOD);
    free(conv_fb);

    /* B(k,n) = sum_{a=0..k} C(k,a) * A^a * S[k-a] */
    i64 *powA = malloc(sizeof(i64) * (k + 1));
    powA[0] = 1;
    for (i64 i = 1; i <= k; i++)
        powA[i] = (i64)((i128)powA[i-1] * A % MOD);

    i64 fk = fact[k];
    i64 ans = 0;
    for (i64 a_cnt = 0; a_cnt <= k; a_cnt++) {
        i64 comb = (i64)((i128)fk * inv_fact[a_cnt] % MOD);
        comb = (i64)((i128)comb * inv_fact[k - a_cnt] % MOD);
        i64 term = (i64)((i128)comb * powA[a_cnt] % MOD);
        term = (i64)((i128)term * S[k - a_cnt] % MOD);
        ans = (ans + term) % MOD;
    }

    free(fact); free(inv_fact); free(fact16); free(inv_fact16);
    free(f); free(fprime); free(b); free(S); free(powA);
    ctx_free(&ctx);

    return (long long)ans;
}
