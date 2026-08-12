/*
 * Project Euler 636 - Restricted Factorisations
 * Native C helper ported from the Python reference solver.
 *
 * Counts representations of n! of the form:
 *   a^1 * b1^2 * b2^2 * c1^3 * c2^3 * c3^3 * d1^4 * d2^4 * d3^4 * d4^4
 * with pairwise distinct bases, order within equal-exponent groups
 * not distinguished.  Uses inclusion-exclusion over set partitions of
 * 10 slots, a coin-change DP for small exponents, and a polynomial
 * quotient-ring recurrence for large exponents.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t  i64;
typedef __int128 i128;

#define MOD       1000000007LL
#define D         30          /* total weight sum */
#define NSLOTS    10
#define MAXBIT    20          /* max exponent < 2^20 for n <= 1e6 */

static const int SLOT_WEIGHTS[NSLOTS] = {1,2,2,3,3,3,4,4,4,4};

/* ---------- modular arithmetic ---------- */

static i64 mod_add(i64 a, i64 b) { i64 r = a + b; if (r >= MOD) r -= MOD; return r; }
static i64 mod_sub(i64 a, i64 b) { i64 r = a - b; if (r < 0) r += MOD; return r; }
static i64 mod_mul(i64 a, i64 b) { return (i64)((i128)a * b % MOD); }

static i64 mod_pow(i64 a, i64 e) {
    i64 r = 1; a %= MOD; if (a < 0) a += MOD;
    while (e > 0) {
        if (e & 1) r = mod_mul(r, a);
        a = mod_mul(a, a);
        e >>= 1;
    }
    return r;
}

/* ---------- set partitions of {0..9} ---------- */

#define MAX_KEYS    60000
#define MAX_KEY_LEN 10

static i64 g_keys[MAX_KEYS][MAX_KEY_LEN];
static int g_key_lens[MAX_KEYS];
static i64 g_coeffs[MAX_KEYS];
static int g_nkeys = 0;

static int g_block_members[10][10];
static int g_block_sizes[10];
static int g_n_blocks = 0;

static int factorial_small(int n) {
    int r = 1;
    for (int i = 2; i <= n; i++) r *= i;
    return r;
}

static void process_partition(void) {
    i64 mu = 1;
    int sums[10];
    int ns = 0;
    for (int b = 0; b < g_n_blocks; b++) {
        int bs = g_block_sizes[b];
        int s = 0;
        for (int m = 0; m < bs; m++)
            s += SLOT_WEIGHTS[g_block_members[b][m]];
        sums[ns++] = s;
        int sign = ((bs - 1) % 2 == 0) ? 1 : -1;
        mu *= (i64)sign * factorial_small(bs - 1);
    }
    /* insertion sort (ns <= 10) */
    for (int i = 1; i < ns; i++) {
        int v = sums[i], j = i - 1;
        while (j >= 0 && sums[j] > v) { sums[j + 1] = sums[j]; j--; }
        sums[j + 1] = v;
    }
    /* find or insert key */
    int found = -1;
    for (int k = 0; k < g_nkeys; k++) {
        if (g_key_lens[k] != ns) continue;
        int match = 1;
        for (int j = 0; j < ns; j++)
            if (g_keys[k][j] != sums[j]) { match = 0; break; }
        if (match) { found = k; break; }
    }
    if (found >= 0) {
        g_coeffs[found] += mu;
    } else {
        int k = g_nkeys++;
        g_key_lens[k] = ns;
        for (int j = 0; j < ns; j++) g_keys[k][j] = sums[j];
        g_coeffs[k] = mu;
    }
}

static void recurse_partition(int i) {
    if (i == NSLOTS) { process_partition(); return; }
    for (int b = 0; b < g_n_blocks; b++) {
        g_block_members[b][g_block_sizes[b]++] = i;
        recurse_partition(i + 1);
        g_block_sizes[b]--;
    }
    g_block_members[g_n_blocks][0] = i;
    g_block_sizes[g_n_blocks] = 1;
    g_n_blocks++;
    recurse_partition(i + 1);
    g_n_blocks--;
}

static void build_partitions(void) {
    g_nkeys = 0;
    g_n_blocks = 0;
    recurse_partition(0);
}

/* ---------- prime sieve ---------- */

static int cmp_int(const void *a, const void *b) {
    int ia = *(const int *)a, ib = *(const int *)b;
    return (ia > ib) - (ia < ib);
}

static int *sieve_primes(int n, int *out_count) {
    char *sieve = (char *)calloc((size_t)n + 1, 1);
    int count = 0;
    for (int i = 2; i <= n; i++) {
        if (!sieve[i]) {
            count++;
            for (long long j = (long long)i * i; j <= n; j += i)
                sieve[j] = 1;
        }
    }
    int *primes = (int *)malloc((size_t)count * sizeof(int));
    int idx = 0;
    for (int i = 2; i <= n; i++)
        if (!sieve[i]) primes[idx++] = i;
    free(sieve);
    *out_count = count;
    return primes;
}

static int factorial_prime_exp(int n, int p) {
    int e = 0;
    long long m = n;
    while (m > 0) { m /= p; e += (int)m; }
    return e;
}

/* ---------- DP: coefficients of prod 1/(1-x^w) up to limit ---------- */

static void coeffs_up_to(const i64 *key, int klen, int limit, i64 *dp) {
    memset(dp, 0, (size_t)(limit + 1) * sizeof(i64));
    dp[0] = 1;
    for (int ki = 0; ki < klen; ki++) {
        int w = (int)key[ki];
        for (int i = w; i <= limit; i++) {
            i64 v = dp[i] + dp[i - w];
            if (v >= MOD) v -= MOD;
            dp[i] = v;
        }
    }
}

/* ---------- polynomial quotient ring (degree D = 30) ---------- */

/* Q(x) = prod(1 - x^w), coefficients not yet reduced mod MOD */
static void poly_Q_from_key(const i64 *key, int klen, i64 *q) {
    memset(q, 0, (size_t)(D + 1) * sizeof(i64));
    q[0] = 1;
    for (int ki = 0; ki < klen; ki++) {
        int w = (int)key[ki];
        for (int i = D - w; i >= 0; i--)
            q[i + w] -= q[i];
    }
}

/* Multiply two degree<D polynomials in the quotient ring.
   r[i] = -q[i+1] mod MOD, i.e. x^D = r[0]*x^{D-1} + ... + r[D-1]. */
static void mul_mod_poly(const i64 *a, const i64 *b, const i64 *r, i64 *out) {
    i64 tmp[2 * D - 1];
    memset(tmp, 0, sizeof(tmp));
    for (int i = 0; i < D; i++) {
        if (!a[i]) continue;
        for (int j = 0; j < D; j++) {
            if (!b[j]) continue;
            tmp[i + j] = mod_add(tmp[i + j], mod_mul(a[i], b[j]));
        }
    }
    for (int k = 2 * D - 2; k >= D; k--) {
        i64 coef = tmp[k];
        if (!coef) continue;
        for (int i = 0; i < D; i++) {
            int idx = k - 1 - i;
            tmp[idx] = mod_add(tmp[idx], mod_mul(coef, r[i]));
        }
    }
    memcpy(out, tmp, (size_t)D * sizeof(i64));
}

static void precompute_x_powers(const i64 *r, int maxbit, i64 *pow_polys) {
    memset(pow_polys, 0, (size_t)maxbit * D * sizeof(i64));
    pow_polys[1] = 1;                       /* x^1 */
    for (int b = 1; b < maxbit; b++)
        mul_mod_poly(pow_polys + (b - 1) * D,
                     pow_polys + (b - 1) * D,
                     r,
                     pow_polys + b * D);
}

static void poly_x_n(const i64 *pow_polys, i64 n, const i64 *r, i64 *out) {
    i64 res[D];
    memset(res, 0, sizeof(res));
    res[0] = 1;                             /* x^0 */
    i64 tmp[D];
    int bit = 0;
    while (n > 0) {
        if (n & 1) {
            memcpy(tmp, res, sizeof(tmp));
            mul_mod_poly(tmp, pow_polys + (size_t)bit * D, r, res);
        }
        n >>= 1;
        bit++;
    }
    memcpy(out, res, sizeof(res));
}

static i64 term_from_poly(const i64 *init, const i64 *poly) {
    i64 s = 0;
    for (int i = 0; i < D; i++)
        s = mod_add(s, mod_mul(poly[i], init[i]));
    return s;
}

/* ---------- main computation F(n!) ---------- */

static i64 compute_F(int n, int cutoff) {
    int nprimes;
    int *primes = sieve_primes(n, &nprimes);

    /* exponent per prime */
    int *exps = (int *)malloc((size_t)nprimes * sizeof(int));
    int max_exp = 0;
    for (int i = 0; i < nprimes; i++) {
        exps[i] = factorial_prime_exp(n, primes[i]);
        if (exps[i] > max_exp) max_exp = exps[i];
    }
    qsort(exps, (size_t)nprimes, sizeof(int), cmp_int);

    /* group into (exponent, count) pairs */
    int *uexp = (int *)malloc((size_t)nprimes * sizeof(int));
    int *ucnt  = (int *)malloc((size_t)nprimes * sizeof(int));
    int nu = 0;
    for (int i = 0; i < nprimes; i++) {
        if (nu > 0 && uexp[nu - 1] == exps[i]) {
            ucnt[nu - 1]++;
        } else {
            uexp[nu] = exps[i];
            ucnt[nu] = 1;
            nu++;
        }
    }

    int use_cutoff = (cutoff < max_exp) ? cutoff : max_exp;
    i64 inv288 = mod_pow(288, MOD - 2);

    int has_exp1 = 0;
    for (int i = 0; i < nu; i++)
        if (uexp[i] == 1) { has_exp1 = 1; break; }

    i64 *dp        = (i64 *)malloc((size_t)(use_cutoff + 1) * sizeof(i64));
    i64 *large_vals= (i64 *)malloc((size_t)nu * sizeof(i64));
    i64  q[D + 1], r[D], init[D], poly_res[D];
    i64  pow_polys[MAXBIT * D];

    i64 total = 0;

    for (int k = 0; k < g_nkeys; k++) {
        int klen = g_key_lens[k];
        const i64 *key = g_keys[k];

        if (has_exp1) {
            int has1 = 0;
            for (int j = 0; j < klen; j++)
                if (key[j] == 1) { has1 = 1; break; }
            if (!has1) continue;
        }

        coeffs_up_to(key, klen, use_cutoff, dp);

        int has_large = 0;
        for (int i = 0; i < nu; i++)
            if (uexp[i] > use_cutoff) { has_large = 1; break; }

        if (has_large) {
            poly_Q_from_key(key, klen, q);
            for (int i = 0; i < D; i++) {
                i64 v = q[i + 1] % MOD;
                if (v < 0) v += MOD;
                r[i] = mod_sub(0, v);
            }
            precompute_x_powers(r, MAXBIT, pow_polys);

            for (int i = 0; i < D; i++)
                init[i] = (i <= use_cutoff) ? dp[i] : 0;

            for (int i = 0; i < nu; i++) {
                if (uexp[i] > use_cutoff) {
                    poly_x_n(pow_polys, uexp[i], r, poly_res);
                    large_vals[i] = term_from_poly(init, poly_res);
                }
            }
        }

        i64 prod = 1;
        for (int i = 0; i < nu; i++) {
            i64 val = (uexp[i] <= use_cutoff) ? dp[uexp[i]] : large_vals[i];
            if (val == 0) { prod = 0; break; }
            prod = mod_mul(prod, mod_pow(val, ucnt[i]));
        }

        if (prod) {
            i64 c = g_coeffs[k] % MOD;
            if (c < 0) c += MOD;
            total = mod_add(total, mod_mul(c, prod));
        }
    }

    total = mod_mul(total, inv288);

    free(large_vals);
    free(dp);
    free(ucnt);
    free(uexp);
    free(exps);
    free(primes);
    return total;
}

long long p636_native(void) {
    build_partitions();
    return (long long)compute_F(1000000, 13000);
}
