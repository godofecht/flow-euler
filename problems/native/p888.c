// Project Euler 888: Coin Game.
// Count losing positions (XOR of Grundy numbers = 0) with m piles, each size [1..N].
// Uses a 16-point ±1 filter over 4 bits of Grundy values.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

static i64 mulmod(i64 a, i64 b, i64 m) {
    return (i64)((i128)a * b % m);
}

static i64 powmod(i64 a, i64 e, i64 m) {
    i64 r = 1 % m; a %= m; if (a < 0) a += m;
    while (e > 0) {
        if (e & 1) r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}

// Extended gcd / modular inverse
static i64 inv_mod(i64 a, i64 mod) {
    a %= mod; if (a < 0) a += mod;
    i64 x0 = 1, x1 = 0, aa = a, mm = mod;
    while (mm) {
        i64 q = aa / mm;
        i64 t = aa - q * mm; aa = mm; mm = t;
        t = x0 - q * x1; x0 = x1; x1 = t;
    }
    if (aa != 1) return -1;
    x0 %= mod; if (x0 < 0) x0 += mod;
    return x0;
}

// Factor n into prime powers
typedef struct { i64 p, pe; } PrimePower;

static int factor_prime_powers(i64 n, PrimePower *out) {
    int count = 0;
    i64 x = n;
    for (i64 p = 2; p * p <= x; p += (p == 2 ? 1 : 2)) {
        if (x % p == 0) {
            i64 pe = 1;
            while (x % p == 0) { x /= p; pe *= p; }
            out[count].p = p;
            out[count].pe = pe;
            count++;
        }
    }
    if (x > 1) {
        out[count].p = x;
        out[count].pe = x;
        count++;
    }
    return count;
}

// Strip p from x: return (x_without_p, v)
static void strip_p(i64 x, i64 p, i64 *x_free, int *v) {
    *x_free = x;
    *v = 0;
    while (*x_free % p == 0) { *x_free /= p; (*v)++; }
}

static const int REMOVE[] = {1, 2, 4, 9};

// Compute Grundy numbers g[0..limit]
static void compute_grundy(int *g, int limit, int split_limit) {
    for (int n = 1; n <= limit; n++) {
        unsigned seen = 0;
        for (int r = 0; r < 4; r++) {
            if (n >= REMOVE[r])
                seen |= 1u << g[n - REMOVE[r]];
        }
        int max_i = (split_limit > 0 && split_limit < n - 1) ? split_limit : n - 1;
        for (int i = 1; i <= max_i; i++)
            seen |= 1u << (g[i] ^ g[n - i]);
        int mex = 0;
        while ((seen >> mex) & 1) mex++;
        g[n] = mex;
    }
}

// Precompute denominator tables for prime power
typedef struct {
    i64 p, pe;
    i64 *den_free, *den_v, *inv_free;
} DenTable;

static void precompute_den_tables(DenTable *dt, int m) {
    i64 p = dt->p, pe = dt->pe;
    dt->den_free = (i64 *)malloc((m + 1) * sizeof(i64));
    dt->den_v = (i64 *)malloc((m + 1) * sizeof(i64));
    dt->inv_free = (i64 *)malloc((m + 1) * sizeof(i64));
    dt->den_free[0] = 1;
    dt->inv_free[0] = 1;
    for (int i = 1; i <= m; i++) {
        i64 xf; int v;
        strip_p(i, p, &xf, &v);
        dt->den_free[i] = xf;
        dt->den_v[i] = v;
        dt->inv_free[i] = inv_mod(xf % pe, pe);
    }
}

// Compute coefficient of x^m in (1-x)^(-a) * (1+x)^(-(N-a)) mod p^e
static i64 coeff_term_mod_prime_power(i64 N, int m, i64 a, i64 p, i64 pe,
                                       const DenTable *dt) {
    i64 A = a, B = N - a;

    i64 *u = (i64 *)calloc(m + 1, sizeof(i64));
    i64 *w = (i64 *)calloc(m + 1, sizeof(i64));

    u[0] = 1 % pe;
    if (A != 0) {
        i64 res = 1 % pe;
        int exp = 0;
        for (int i = 1; i <= m; i++) {
            i64 num = A + i - 1;
            i64 num_free; int v_num;
            strip_p(num, p, &num_free, &v_num);
            exp += v_num;
            exp -= (int)dt->den_v[i];
            res = mulmod(res, num_free % pe, pe);
            res = mulmod(res, dt->inv_free[i], pe);
            u[i] = mulmod(res, powmod(p, exp, pe), pe);
        }
    }

    w[0] = 1 % pe;
    if (B != 0) {
        i64 res = 1 % pe;
        int exp = 0;
        for (int i = 1; i <= m; i++) {
            i64 num = B + i - 1;
            i64 num_free; int v_num;
            strip_p(num, p, &num_free, &v_num);
            exp += v_num;
            exp -= (int)dt->den_v[i];
            res = mulmod(res, num_free % pe, pe);
            res = mulmod(res, dt->inv_free[i], pe);
            i64 val = mulmod(res, powmod(p, exp, pe), pe);
            if (i & 1) val = (pe - val) % pe;
            w[i] = val;
        }
    }

    i64 out = 0;
    for (int i = 0; i <= m; i++)
        out = (out + mulmod(u[i], w[m - i], pe)) % pe;

    free(u);
    free(w);
    return out;
}

// Count Grundy values in [1..N] using periodicity
static void grundy_counts_up_to(i64 N, const int *g, int pre, int period,
                                 i64 *counts, int max_g) {
    memset(counts, 0, max_g * sizeof(i64));
    if (N <= 0) return;

    if (N < pre) {
        for (i64 n = 1; n <= N; n++)
            counts[g[n]]++;
        return;
    }

    for (int n = 1; n < pre; n++)
        counts[g[n]]++;

    i64 per_counts[16] = {0};
    for (int n = pre; n < pre + period; n++)
        per_counts[g[n]]++;

    i64 total_period_terms = N - pre + 1;
    i64 q = total_period_terms / period;
    i64 r = total_period_terms % period;

    for (int k = 0; k < max_g; k++)
        counts[k] += per_counts[k] * q;

    for (i64 n = pre; n < pre + r; n++)
        counts[g[n]]++;
}

static int popcount(int x) { return __builtin_popcount(x); }

static i64 S_mod(i64 N, int m, i64 mod, const int *g, int pre, int period) {
    int max_g = 16;
    i64 counts[16];
    grundy_counts_up_to(N, g, pre, period, counts, max_g);

    PrimePower pp[20];
    int npp = factor_prime_powers(mod, pp);

    i64 mods[20], Ms[20], inv_Ms[20];
    for (int i = 0; i < npp; i++) {
        mods[i] = pp[i].pe;
        Ms[i] = mod / pp[i].pe;
        inv_Ms[i] = inv_mod(Ms[i] % pp[i].pe, pp[i].pe);
    }

    DenTable den_tables[20];
    for (int i = 0; i < npp; i++) {
        den_tables[i].p = pp[i].p;
        den_tables[i].pe = pp[i].pe;
        precompute_den_tables(&den_tables[i], m);
    }

    i64 inv16 = inv_mod(16, mod);
    i64 total = 0;

    for (int s = 0; s < 16; s++) {
        i64 a = 0;
        for (int gv = 0; gv < max_g; gv++) {
            if (counts[gv] && (popcount(gv & s) & 1) == 0)
                a += counts[gv];
        }

        // CRT combine
        i64 x = 0;
        for (int i = 0; i < npp; i++) {
            i64 resid = coeff_term_mod_prime_power(N, m, a,
                den_tables[i].p, den_tables[i].pe, &den_tables[i]);
            x = (x + mulmod(mulmod(resid % mods[i], Ms[i], mod), inv_Ms[i], mod)) % mod;
        }

        total = (total + x) % mod;
    }

    for (int i = 0; i < npp; i++) {
        free(den_tables[i].den_free);
        free(den_tables[i].den_v);
        free(den_tables[i].inv_free);
    }

    return mulmod(total, inv16, mod);
}

long long p888_native(void) {
    i64 N = 12491249;
    int m = 1249;
    i64 MOD = 912491249;

    int PRE = 322;
    int PERIOD = 11060;
    int SPLIT_LIMIT = 600;

    int precomp_limit = PRE + 2 * PERIOD;
    int *g_fast = (int *)malloc((precomp_limit + 1) * sizeof(int));
    g_fast[0] = 0;
    compute_grundy(g_fast, precomp_limit, SPLIT_LIMIT);

    i64 ans = S_mod(N, m, MOD, g_fast, PRE, PERIOD);

    free(g_fast);
    return ans;
}
