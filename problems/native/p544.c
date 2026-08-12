/* Project Euler 544 - Chromatic Conundrum
 * F(r,c,n) = proper colourings of r x c grid with at most n colours.
 * S(r,c,n) = sum_{k=1..n} F(r,c,k).
 * Compute S(9, 10, 1112131415) mod 1e9+7.
 *
 * Frontier DP for chromatic polynomial, then Lagrange interpolation.
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MOD 1000000007LL
#define R 9
#define C 10
#define MAXSTATES 65536
#define MAXPOLY 128

typedef int64_t i64;
typedef __int128 i128;

static i64 mod_pow(i64 base, i64 exp, i64 mod) {
    i64 r = 1 % mod;
    base %= mod;
    if (base < 0) base += mod;
    while (exp > 0) {
        if (exp & 1) r = (i128)r * base % mod;
        base = (i128)base * base % mod;
        exp >>= 1;
    }
    return r;
}

static i64 mod_inv(i64 a, i64 mod) {
    return mod_pow(a, mod - 2, mod);
}

/* Polynomial in q (power basis), coefficients mod MOD. */
typedef struct {
    int len;
    i64 coeff[MAXPOLY];
} Poly;

static void poly_zero(Poly *p) {
    p->len = 0;
    memset(p->coeff, 0, sizeof(p->coeff));
}

static void poly_add(Poly *tgt, const Poly *src) {
    if (src->len > tgt->len) {
        for (int i = tgt->len; i < src->len; i++) tgt->coeff[i] = 0;
        tgt->len = src->len;
    }
    for (int i = 0; i < src->len; i++) {
        tgt->coeff[i] += src->coeff[i];
        if (tgt->coeff[i] >= MOD) tgt->coeff[i] -= MOD;
    }
}

/* tgt += (q - m) * poly */
static void poly_add_q_minus_m(Poly *tgt, const Poly *poly, i64 m) {
    int need = poly->len + 1;
    if (need > tgt->len) {
        for (int i = tgt->len; i < need; i++) tgt->coeff[i] = 0;
        tgt->len = need;
    }
    i64 mm = m % MOD;
    for (int j = 0; j < poly->len; j++) {
        i64 v = poly->coeff[j];
        /* coeff[j] -= m * v */
        i64 sub = (i128)mm * v % MOD;
        tgt->coeff[j] = (tgt->coeff[j] - sub + MOD) % MOD;
        /* coeff[j+1] += v */
        tgt->coeff[j + 1] += v;
        if (tgt->coeff[j + 1] >= MOD) tgt->coeff[j + 1] -= MOD;
    }
}

static i64 poly_eval(const Poly *p, i64 x) {
    x %= MOD;
    if (x < 0) x += MOD;
    i64 res = 0, pw = 1;
    for (int i = 0; i < p->len; i++) {
        res = (res + (i128)p->coeff[i] * pw) % MOD;
        pw = (i128)pw * x % MOD;
    }
    return res;
}

/* State: tuple of R labels, -1 = EMPTY. Canonical form. */
typedef struct {
    int8_t labels[R];
} State;

static int state_eq(const State *a, const State *b) {
    return memcmp(a->labels, b->labels, R) == 0;
}

static int state_hash(const State *s) {
    unsigned int h = 0;
    for (int i = 0; i < R; i++) {
        h = h * 31 + (unsigned int)(s->labels[i] + 2);
    }
    return (int)(h & (MAXSTATES - 1));
}

/* Canonical relabel by first appearance. */
static State canon(const State *s) {
    State out;
    int8_t mp[R + 2];
    memset(mp, -1, sizeof(mp));
    int nxt = 0;
    for (int i = 0; i < R; i++) {
        if (s->labels[i] == -1) {
            out.labels[i] = -1;
        } else {
            int8_t y = mp[s->labels[i] + 1];
            if (y == -1) {
                y = (int8_t)nxt;
                mp[s->labels[i] + 1] = y;
                nxt++;
            }
            out.labels[i] = y;
        }
    }
    return out;
}

/* DP map: state -> Poly. Use hash table. */
typedef struct {
    State key;
    Poly val;
    int used;
} Entry;

static void dp_init(Entry *table) {
    memset(table, 0, sizeof(Entry) * MAXSTATES);
}

static Poly *dp_get_or_create(Entry *table, const State *key) {
    int h = state_hash(key);
    for (int i = 0; i < MAXSTATES; i++) {
        int idx = (h + i) & (MAXSTATES - 1);
        if (!table[idx].used) {
            table[idx].key = *key;
            table[idx].used = 1;
            poly_zero(&table[idx].val);
            return &table[idx].val;
        }
        if (state_eq(&table[idx].key, key)) {
            return &table[idx].val;
        }
    }
    /* Table full - should not happen */
    fprintf(stderr, "DP table full!\n");
    exit(1);
    return NULL;
}

static void dp_copy(Entry *dst, Entry *src) {
    memcpy(dst, src, sizeof(Entry) * MAXSTATES);
}

static void dp_clear(Entry *table) {
    memset(table, 0, sizeof(Entry) * MAXSTATES);
}

/* Compute chromatic polynomial coefficients for r x c grid. */
static void chromatic_poly(int r, int c, Poly *out) {
    static Entry dp[MAXSTATES], ndp[MAXSTATES];
    dp_init(dp);

    State start;
    memset(&start, -1, sizeof(State));
    Poly *p = dp_get_or_create(dp, &start);
    p->len = 1;
    p->coeff[0] = 1;

    for (int col = 0; col < c; col++) {
        for (int row = 0; row < r; row++) {
            int i = row;
            int has_up = row > 0;
            int has_left = col > 0;
            dp_clear(ndp);

            for (int idx = 0; idx < MAXSTATES; idx++) {
                if (!dp[idx].used) continue;
                State st = dp[idx].key;
                Poly *poly = &dp[idx].val;

                /* Compute m = number of distinct labels */
                int mx = -1;
                for (int j = 0; j < r; j++) {
                    if (st.labels[j] > mx) mx = st.labels[j];
                }
                int m = mx + 1;

                int8_t up = has_up ? st.labels[i - 1] : -2;
                int8_t left = has_left ? st.labels[i] : -2;

                /* Try each existing label */
                for (int lab = 0; lab < m; lab++) {
                    if (lab == up || lab == left) continue;
                    State ns = st;
                    ns.labels[i] = (int8_t)lab;
                    ns = canon(&ns);
                    Poly *tgt = dp_get_or_create(ndp, &ns);
                    if (tgt) poly_add(tgt, poly);
                }

                /* New label (q - m) choices */
                State ns = st;
                ns.labels[i] = (int8_t)m;
                ns = canon(&ns);
                Poly *tgt = dp_get_or_create(ndp, &ns);
                if (tgt) poly_add_q_minus_m(tgt, poly, m);
            }
            dp_copy(dp, ndp);
        }
    }

    /* Forget frontier */
    for (int i = 0; i < r; i++) {
        dp_clear(ndp);
        for (int idx = 0; idx < MAXSTATES; idx++) {
            if (!dp[idx].used) continue;
            State st = dp[idx].key;
            st.labels[i] = -1;
            State ns = canon(&st);
            Poly *tgt = dp_get_or_create(ndp, &ns);
            if (tgt) poly_add(tgt, &dp[idx].val);
        }
        dp_copy(dp, ndp);
    }

    /* Find start state */
    *out = *dp_get_or_create(dp, &start);
}

/* Lagrange interpolation for x_i = 0..n */
static i64 lagrange_eval(const i64 *values, int n, i64 x) {
    if (x >= 0 && x <= n) return values[x] % MOD;
    x %= MOD;
    if (x < 0) x += MOD;

    i64 fac[n + 1];
    fac[0] = 1;
    for (int i = 1; i <= n; i++) fac[i] = (i128)fac[i - 1] * i % MOD;

    i64 invfac[n + 1];
    invfac[n] = mod_inv(fac[n], MOD);
    for (int i = n; i > 0; i--) invfac[i - 1] = (i128)invfac[i] * i % MOD;

    i64 pre[n + 2], suf[n + 2];
    pre[0] = 1;
    for (int i = 0; i <= n; i++) {
        i64 term = (x - i) % MOD;
        if (term < 0) term += MOD;
        pre[i + 1] = (i128)pre[i] * term % MOD;
    }
    suf[n + 1] = 1;
    for (int i = n; i >= 0; i--) {
        i64 term = (x - i) % MOD;
        if (term < 0) term += MOD;
        suf[i] = (i128)suf[i + 1] * term % MOD;
    }

    i64 ans = 0;
    for (int i = 0; i <= n; i++) {
        i64 num = (i128)pre[i] * suf[i + 1] % MOD;
        i64 den = (i128)invfac[i] * invfac[n - i] % MOD;
        if ((n - i) & 1) den = (MOD - den) % MOD;
        i64 term = (i128)values[i] * num % MOD;
        term = (i128)term * den % MOD;
        ans = (ans + term) % MOD;
    }
    return ans;
}

long long p544_native(void) {
    Poly p;
    chromatic_poly(R, C, &p);

    int deg = R * C;
    int m = deg + 1;
    i64 prefix[m + 2];
    prefix[0] = 0;
    i64 acc = 0;
    for (int k = 1; k <= m + 1; k++) {
        acc = (acc + poly_eval(&p, k)) % MOD;
        prefix[k] = acc;
    }

    i64 n = 1112131415LL;
    i64 ans = lagrange_eval(prefix, m, n);
    return (long long)ans;
}
