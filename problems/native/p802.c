/* Project Euler 802: Iterated Composition.
 *
 * The map f(x,y) = (x^2 - x - y^2, 2xy - y + pi) corresponds to the complex
 * polynomial F(z) = z^2 - z + i*pi.  The requested sum is always an integer and
 * can be computed purely arithmetically.
 *
 * Let A(d) be the sum of x-coordinates of points with period dividing d.
 *   A(1) = 2,  A(d) = 2^(d-1) for d >= 2.
 * Let S(d) be the sum for exact period d.  Then A(n) = sum_{d|n} S(d), so by
 * Mobius inversion:
 *   P(n) = sum_{d<=n} A(d) * M(floor(n/d))
 * where M is the Mertens function.  Evaluate with floor-division grouping.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef int64_t i64;

#define MOD 1020340567LL
#define N_DEFAULT 10000000LL

static i64 mod_pow(i64 a, i64 e, i64 mod) {
    i64 r = 1 % mod; a %= mod;
    while (e > 0) {
        if (e & 1) r = (i64)((__int128)r * a % mod);
        a = (i64)((__int128)a * a % mod);
        e >>= 1;
    }
    return r;
}

/* Build sorted distinct values of floor(n / k) for k=1..n.
 * Returns the count via *out_count; caller frees the buffer. */
static i64 *build_floor_div_queries(i64 n, i64 *out_count) {
    i64 cap = 64, cnt = 0;
    i64 *qs = (i64 *)malloc(cap * sizeof(i64));
    i64 k = 1;
    while (k <= n) {
        i64 q = n / k;
        if (cnt >= cap) { cap *= 2; qs = (i64 *)realloc(qs, cap * sizeof(i64)); }
        qs[cnt++] = q;
        k = n / q + 1;
    }
    /* qs is already strictly decreasing; reverse to ascending. */
    for (i64 i = 0, j = cnt - 1; i < j; i++, j--) {
        i64 t = qs[i]; qs[i] = qs[j]; qs[j] = t;
    }
    *out_count = cnt;
    return qs;
}

/* Compute Mertens M(t) = sum_{m<=t} mu(m) for the selected query points.
 * Linear sieve for mu up to n.  Returns a parallel array of Mertens values
 * matching the ascending query array. */
static i64 *mertens_at_points(i64 n, const i64 *points, i64 npoints) {
    if (npoints == 0) return NULL;
    i64 *res = (i64 *)calloc(npoints, sizeof(i64));

    signed char *mu = (signed char *)calloc((size_t)n + 1, 1);
    i64 *lp = (i64 *)calloc((size_t)n + 1, sizeof(i64));
    i64 *primes = (i64 *)malloc((size_t)n * sizeof(i64));
    i64 nprimes = 0;

    mu[1] = 1;
    i64 mertens = 1;

    i64 idx = 0;
    if (points[0] == 1) { res[0] = 1; idx = 1; }

    for (i64 i = 2; i <= n; i++) {
        if (lp[i] == 0) {
            lp[i] = i;
            primes[nprimes++] = i;
            mu[i] = -1;
        }
        i64 li = lp[i];
        signed char mui = mu[i];
        for (i64 pi = 0; pi < nprimes; pi++) {
            i64 p = primes[pi];
            if (p > li) break;
            i64 ip = i * p;
            if (ip > n) break;
            lp[ip] = p;
            if (p == li) { mu[ip] = 0; break; }
            else { mu[ip] = (signed char)(-mui); }
        }
        mertens += mu[i];
        if (idx < npoints && i == points[idx]) {
            res[idx] = mertens;
            idx++;
        }
    }

    free(mu);
    free(lp);
    free(primes);
    return res;
}

/* Sum_{d=l..r} A(d) mod MOD, where A(1)=2, A(d)=2^(d-1) for d>=2.
 *   sum_{d=1..r} A(d) = 2^r
 *   for l>=2: sum_{d=l..r} A(d) = 2^r - 2^(l-1)
 */
static i64 sum_A(i64 l, i64 r) {
    if (l == 1) return mod_pow(2, r, MOD);
    i64 v = (mod_pow(2, r, MOD) - mod_pow(2, l - 1, MOD)) % MOD;
    if (v < 0) v += MOD;
    return v;
}

static i64 P_mod(i64 n) {
    if (n <= 0) return 0;
    i64 nq;
    i64 *qs = build_floor_div_queries(n, &nq);
    i64 *mertens = mertens_at_points(n, qs, nq);

    i64 ans = 0;
    i64 l = 1;
    while (l <= n) {
        i64 q = n / l;
        i64 r = n / q;
        /* find index of q in ascending qs via binary search */
        i64 lo = 0, hi = nq - 1, found = -1;
        while (lo <= hi) {
            i64 mid = (lo + hi) / 2;
            if (qs[mid] == q) { found = mid; break; }
            else if (qs[mid] < q) lo = mid + 1;
            else hi = mid - 1;
        }
        i64 m = mertens[found]; /* can be negative */
        i64 s = sum_A(l, r);
        i64 term = (i64)((__int128)s * m % MOD);
        ans = (ans + term) % MOD;
        if (ans < 0) ans += MOD;
        l = r + 1;
    }

    free(qs);
    free(mertens);
    return ans;
}

long long p802_native(void) {
    return (long long)P_mod(N_DEFAULT);
}
