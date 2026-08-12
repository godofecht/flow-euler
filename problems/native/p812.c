/* Project Euler 812: Dynamical Polynomials
   S(10000) mod 998244353 via NTT-based polynomial exp/ln.
   Port of the Python reference solver. */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef long long ll;
typedef __int128 i128;

#define MOD 998244353LL
#define PRIMITIVE_ROOT 3LL

static ll powmod(ll base, ll exp, ll mod) {
    ll r = 1 % mod;
    base %= mod;
    if (base < 0) base += mod;
    while (exp) {
        if (exp & 1) r = r * base % mod;
        base = base * base % mod;
        exp >>= 1;
    }
    return r;
}

static ll inv_mod(ll a) { return powmod(a, MOD - 2, MOD); }

/* Sieve primes <= n */
static int *sieve_primes(int n, int *out_cnt) {
    if (n < 2) { *out_cnt = 0; return NULL; }
    unsigned char *bs = (unsigned char*)malloc(n + 1);
    memset(bs, 1, n + 1);
    bs[0] = bs[1] = 0;
    for (int p = 2; (ll)p * p <= n; p++) {
        if (bs[p]) {
            for (int j = p * p; j <= n; j += p) bs[j] = 0;
        }
    }
    int cnt = 0;
    for (int i = 2; i <= n; i++) if (bs[i]) cnt++;
    int *primes = (int*)malloc(cnt * sizeof(int));
    int idx = 0;
    for (int i = 2; i <= n; i++) if (bs[i]) primes[idx++] = i;
    free(bs);
    *out_cnt = cnt;
    return primes;
}

/* In-place iterative NTT */
static void ntt(ll *a, int n, int invert) {
    int j = 0;
    for (int i = 1; i < n; i++) {
        int bit = n >> 1;
        while (j & bit) { j ^= bit; bit >>= 1; }
        j ^= bit;
        if (i < j) { ll tmp = a[i]; a[i] = a[j]; a[j] = tmp; }
    }

    for (int len = 2; len <= n; len <<= 1) {
        ll wlen = powmod(PRIMITIVE_ROOT, (MOD - 1) / len, MOD);
        if (invert) wlen = powmod(wlen, MOD - 2, MOD);
        int half = len >> 1;
        for (int i = 0; i < n; i += len) {
            ll w = 1;
            for (int j2 = 0; j2 < half; j2++) {
                ll u = a[i + j2];
                ll v = a[i + j2 + half] * w % MOD;
                a[i + j2] = (u + v) % MOD;
                a[i + j2 + half] = ((u - v) % MOD + MOD) % MOD;
                w = w * wlen % MOD;
            }
        }
    }

    if (invert) {
        ll inv_n = powmod(n, MOD - 2, MOD);
        for (int i = 0; i < n; i++) a[i] = a[i] * inv_n % MOD;
    }
}

/* Polynomial multiply mod MOD */
static ll *polymul(const ll *a, int la, const ll *b, int lb, int *out_len) {
    if (la == 0 || lb == 0) { *out_len = 0; return NULL; }
    int n = la + lb - 1;

    /* Small-size fallback */
    if (la <= 32 || lb <= 32) {
        ll *res = (ll*)calloc(n, sizeof(ll));
        for (int i = 0; i < la; i++) {
            if (a[i] == 0) continue;
            for (int j = 0; j < lb; j++) {
                res[i + j] = (res[i + j] + a[i] * b[j]) % MOD;
            }
        }
        *out_len = n;
        return res;
    }

    int size = 1;
    while (size < n) size <<= 1;
    ll *fa = (ll*)calloc(size, sizeof(ll));
    ll *fb = (ll*)calloc(size, sizeof(ll));
    memcpy(fa, a, la * sizeof(ll));
    memcpy(fb, b, lb * sizeof(ll));
    ntt(fa, size, 0);
    ntt(fb, size, 0);
    for (int i = 0; i < size; i++) fa[i] = fa[i] * fb[i] % MOD;
    ntt(fa, size, 1);
    ll *res = (ll*)malloc(n * sizeof(ll));
    memcpy(res, fa, n * sizeof(ll));
    free(fa); free(fb);
    *out_len = n;
    return res;
}

static ll *poly_der(const ll *a, int la, int *out_len) {
    if (la <= 1) { *out_len = 0; return NULL; }
    *out_len = la - 1;
    ll *res = (ll*)malloc(*out_len * sizeof(ll));
    for (int i = 1; i < la; i++) res[i - 1] = (ll)i * a[i] % MOD;
    return res;
}

static ll *poly_int(const ll *a, int la, const ll *invs, int *out_len) {
    *out_len = la + 1;
    ll *res = (ll*)calloc(*out_len, sizeof(ll));
    for (int i = 0; i < la; i++) res[i + 1] = a[i] * invs[i + 1] % MOD;
    return res;
}

/* Inverse series of a mod x^n. Requires a[0] != 0. */
static ll *poly_inv(const ll *a, int la, int n, int *out_len) {
    *out_len = n;
    ll *res = (ll*)calloc(n, sizeof(ll));
    res[0] = powmod(a[0], MOD - 2, MOD);
    int m = 1;
    while (m < n) {
        int m2 = m * 2;
        if (m2 > n) m2 = n;
        int t_len;
        ll *t = polymul(a, la > m2 ? m2 : la, res, m, &t_len);
        /* t has length min(la, m2) + m - 1, but we only need first m2 */
        /* t = 2 - t (mod) */
        int tl = t_len > m2 ? m2 : t_len;
        ll *t2 = (ll*)calloc(m2, sizeof(ll));
        for (int i = 0; i < tl; i++) t2[i] = t[i];
        free(t);
        t2[0] = (2 - t2[0]) % MOD;
        if (t2[0] < 0) t2[0] += MOD;
        for (int i = 1; i < m2; i++) { t2[i] = (-t2[i]) % MOD; if (t2[i] < 0) t2[i] += MOD; }
        int r_len;
        ll *new_res = polymul(res, m, t2, m2, &r_len);
        int rl = r_len > m2 ? m2 : r_len;
        memcpy(res, new_res, rl * sizeof(ll));
        free(new_res); free(t2);
        m = m2;
    }
    return res;
}

/* ln(a) mod x^n. Requires a[0] == 1. */
static ll *poly_ln(const ll *a, int la, int n, const ll *invs, int *out_len) {
    int der_len;
    ll *der = poly_der(a, la, &der_len);
    int inv_len;
    ll *inv_a = poly_inv(a, la, n, &inv_len);
    int q_len;
    ll *q = polymul(der, der_len, inv_a, inv_len, &q_len);
    int ql = n - 1 > q_len ? q_len : (n - 1);
    if (ql < 0) ql = 0;
    int int_len;
    ll *result = poly_int(q, ql, invs, &int_len);
    int rl = n > int_len ? int_len : n;
    *out_len = rl;
    free(der); free(inv_a); free(q);
    /* ensure length n */
    ll *res = (ll*)calloc(n, sizeof(ll));
    memcpy(res, result, rl * sizeof(ll));
    free(result);
    return res;
}

/* exp(f) mod x^n. Requires f[0] == 0. */
static ll *poly_exp(const ll *f, int lf, int n, const ll *invs, int *out_len) {
    *out_len = n;
    ll *g = (ll*)calloc(n, sizeof(ll));
    g[0] = 1;
    int m = 1;
    while (m < n) {
        int m2 = m * 2;
        if (m2 > n) m2 = n;
        /* g_pad = g padded to m2 */
        ll *g_pad = (ll*)calloc(m2, sizeof(ll));
        memcpy(g_pad, g, m * sizeof(ll));
        int ln_len;
        ll *ln_g = poly_ln(g_pad, m2, m2, invs, &ln_len);
        ll *diff = (ll*)calloc(m2, sizeof(ll));
        for (int i = 0; i < m2; i++) {
            ll fi = (i < lf) ? f[i] : 0;
            diff[i] = ((fi - ln_g[i]) % MOD + MOD) % MOD;
        }
        diff[0] = (diff[0] + 1) % MOD;
        int new_len;
        ll *new_g = polymul(g, m, diff, m2, &new_len);
        int nl = new_len > m2 ? m2 : new_len;
        memcpy(g, new_g, nl * sizeof(ll));
        free(g_pad); free(ln_g); free(diff); free(new_g);
        m = m2;
    }
    return g;
}

/* Special component generating function */
static ll *special_component(int N) {
    ll inv2 = (MOD + 1) / 2;
    ll *v1 = (ll*)calloc(N + 1, sizeof(ll));
    v1[0] = 1;
    for (int p = 2; p <= N; p <<= 1) {
        for (int d = p; d <= N; d++) {
            v1[d] = (v1[d] + v1[d - p]) % MOD;
        }
    }

    ll *vminus = (ll*)calloc(N + 1, sizeof(ll));
    vminus[0] = 1;
    for (int p = 2; p <= N; p <<= 1) {
        ll *new_arr = (ll*)calloc(N + 1, sizeof(ll));
        for (int d = 0; d <= N; d++) {
            ll val = vminus[d];
            if (d >= p) {
                val = ((val - new_arr[d - p]) % MOD + MOD) % MOD;
            }
            new_arr[d] = val;
        }
        free(vminus);
        vminus = new_arr;
    }

    /* P(x) = 1/2 * ((1+x)V1 + (1-x)Vminus) */
    ll *pser = (ll*)calloc(N + 1, sizeof(ll));
    for (int d = 0; d <= N; d++) {
        ll t1 = (v1[d] + (d > 0 ? v1[d - 1] : 0)) % MOD;
        ll t2 = ((vminus[d] - (d > 0 ? vminus[d - 1] : 0)) % MOD + MOD) % MOD;
        pser[d] = ((t1 + t2) % MOD) * inv2 % MOD;
    }

    /* Multiply by 1/(1-x): prefix sums */
    ll *a = (ll*)calloc(N + 1, sizeof(ll));
    ll s = 0;
    for (int d = 0; d <= N; d++) {
        s = (s + pser[d]) % MOD;
        a[d] = s;
    }

    /* Multiply by 1/(1-x^2) */
    ll *b = (ll*)calloc(N + 1, sizeof(ll));
    for (int d = 0; d <= N; d++) {
        ll val = a[d];
        if (d >= 2) val = (val + b[d - 2]) % MOD;
        b[d] = val;
    }

    free(v1); free(vminus); free(pser); free(a);
    return b;
}

/* Add component multiplicities */
static void add_component_multiplicities(ll *c, int N) {
    int limit_phi = 2 * N;
    int nprimes;
    int *primes = sieve_primes(limit_phi + 1, &nprimes);

    /* filter odd primes */
    int *odd_primes = (int*)malloc(nprimes * sizeof(int));
    int nop = 0;
    for (int i = 0; i < nprimes; i++) {
        if (primes[i] & 1) odd_primes[nop++] = primes[i];
    }

    /* DFS using explicit stack to avoid deep recursion */
    /* Stack frame: (start_idx, n_val, phi_val) */
    typedef struct { int start_idx; long long n_val; long long phi_val; } Frame;
    int stack_cap = 200000;
    Frame *stack = (Frame*)malloc(stack_cap * sizeof(Frame));
    int sp = 0;
    stack[sp++] = (Frame){0, 1, 1};

    while (sp > 0) {
        Frame f = stack[--sp];

        /* Record this odd m0 (exclude 1) */
        if (f.n_val > 1) {
            long long ph = f.phi_val;
            long long deg0 = ph / 2;
            long long w = 0;
            int k = 0;
            while (1) {
                long long deg;
                if (k == 0 || k == 1) deg = deg0;
                else deg = ph * (1LL << (k - 2));
                w += deg;
                if (w > N) break;
                c[w] += 1;
                k++;
            }
        }

        /* Extend by adding new prime powers */
        for (int i = f.start_idx; i < nop; i++) {
            int p = odd_primes[i];
            if (f.phi_val * (p - 1) > limit_phi) break;

            /* exponent 1 */
            if (sp >= stack_cap) { stack_cap *= 2; stack = (Frame*)realloc(stack, stack_cap * sizeof(Frame)); }
            stack[sp++] = (Frame){i + 1, f.n_val * p, f.phi_val * (p - 1)};

            /* exponent >= 2 */
            long long n_e = (long long)f.n_val * p * p;
            long long phi_e = (long long)f.phi_val * (p - 1) * p;
            while (phi_e <= limit_phi) {
                if (sp >= stack_cap) { stack_cap *= 2; stack = (Frame*)realloc(stack, stack_cap * sizeof(Frame)); }
                stack[sp++] = (Frame){i + 1, n_e, phi_e};
                n_e *= p;
                phi_e *= p;
            }
        }
    }

    free(stack);
    free(primes); free(odd_primes);
}

static ll solve(int N) {
    ll *c = (ll*)calloc(N + 1, sizeof(ll));
    add_component_multiplicities(c, N);

    /* Build g(x) = log(F)(x) where F(x) = Prod (1 - x^a)^(-c[a]) */
    ll *g = (ll*)calloc(N + 1, sizeof(ll));
    for (int d = 1; d <= N; d++) {
        ll cd = c[d];
        if (cd) {
            ll add = (ll)d * (cd % MOD) % MOD;
            for (int k = d; k <= N; k += d) {
                g[k] = (g[k] + add) % MOD;
            }
        }
    }

    /* Inverses */
    ll *invs = (ll*)calloc(N + 2, sizeof(ll));
    for (int i = 1; i <= N + 1; i++) invs[i] = inv_mod(i);

    /* h[k] = g[k] * invs[k] */
    ll *h = (ll*)calloc(N + 1, sizeof(ll));
    for (int k = 1; k <= N; k++) h[k] = g[k] * invs[k] % MOD;

    int colored_len;
    ll *colored = poly_exp(h, N + 1, N + 1, invs, &colored_len);

    ll *special = special_component(N);

    /* Total = convolution of colored and special */
    int total_len;
    ll *total = polymul(colored, N + 1, special, N + 1, &total_len);

    ll result = total[N] % MOD;

    free(c); free(g); free(invs); free(h); free(colored); free(special); free(total);
    return result;
}

long long p812_native(void) {
    return solve(10000);
}
