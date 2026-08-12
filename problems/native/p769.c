/*
 * Project Euler 769 - Count primitive representations of squares by
 * f(x,y) = x^2 + 5xy + 3y^2, i.e. f(x,y) = z^2 with z <= N, x,y>0, gcd(x,y)=1.
 *
 * Two branches:
 *   1. p > 0, q > sqrt(3)*p, z = q^2 + 5pq + 3p^2
 *   2. p = -a < 0, sqrt(3)*a < q < 2.5a, z = -(q^2 - 5aq + 3a^2)
 *
 * For each p (or a), count coprime q in range via inclusion-exclusion over
 * distinct prime factors.  Non-primitive cases (q ≡ 4p mod 13) are subtracted.
 *
 * Ported from the Python reference solver.
 */

#include <stdint.h>
#include <stdlib.h>
#include <math.h>

typedef int64_t i64;

static const double SQRT3 = 1.7320508075688772;

/* ---- integer square root for i64 ---- */
static i64 isqrt64(i64 n) {
    if (n <= 0) return 0;
    i64 x = (i64)sqrt((double)n);
    while (x > 0 && x * x > n) x--;
    while ((x + 1) * (x + 1) <= n) x++;
    return x;
}

/* ---- inverses mod 13 ---- */
static int inv_mod13[13];

static void init_inv_mod13(void) {
    for (int a = 1; a < 13; a++)
        for (int x = 1; x < 13; x++)
            if ((a * x) % 13 == 1) { inv_mod13[a] = x; break; }
}

/* ---- linear sieve for smallest prime factor ---- */
static int *spf_arr = NULL;

static void build_spf(int n) {
    spf_arr = calloc(n + 1, sizeof(int));
    int *primes = malloc((size_t)n * sizeof(int));
    int pc = 0;
    for (int i = 2; i <= n; i++) {
        if (!spf_arr[i]) { spf_arr[i] = i; primes[pc++] = i; }
        for (int j = 0; j < pc; j++) {
            long v = (long)primes[j] * i;
            if (v > n) break;
            spf_arr[v] = primes[j];
            if (primes[j] == spf_arr[i]) break;
        }
    }
    free(primes);
}

/* ---- distinct prime factors of n (n <= spf table size) ---- */
static int distinct_prime_factors(i64 n, int *out) {
    int cnt = 0;
    while (n > 1) {
        int p = spf_arr[n];
        out[cnt++] = p;
        while (n % p == 0) n /= p;
    }
    return cnt;
}

/* ---- generate squarefree divisors and Mobius values ---- */
static int gen_divisors_mu(const int *primes, int np, i64 *ds, int *mus) {
    ds[0] = 1; mus[0] = 1;
    int len = 1;
    for (int i = 0; i < np; i++) {
        int p = primes[i];
        for (int j = 0; j < len; j++) {
            ds[j + len] = ds[j] * p;
            mus[j + len] = -mus[j];
        }
        len *= 2;
    }
    return len;
}

/* ---- count integers in [L, R] with x ≡ rem (mod mod), 0 <= rem < mod ---- */
static i64 count_cong(i64 L, i64 R, i64 mod, i64 rem) {
    if (rem < L)
        rem += ((L - rem + mod - 1) / mod) * mod;
    if (rem > R) return 0;
    return 1 + (R - rem) / mod;
}

/* ---- count coprime q in [L, R] with gcd(q, n) = 1 ---- */
static i64 count_coprime_interval(i64 L, i64 R, const i64 *ds, const int *mus, int nd) {
    i64 Lm = L - 1;
    i64 total = 0;
    for (int i = 0; i < nd; i++)
        total += (i64)mus[i] * (R / ds[i] - Lm / ds[i]);
    return total;
}

/* ---- count coprime q in [L, R] with q ≡ rem13 (mod 13) ---- */
static i64 count_coprime_mod13(const i64 *ds, const int *mus, int nd,
                                i64 L, i64 R, int rem13) {
    i64 total = 0;
    for (int i = 0; i < nd; i++) {
        i64 d = ds[i];
        int inv = inv_mod13[(int)(d % 13)];
        int m0 = (rem13 * inv) % 13;
        i64 rem = d * m0;
        i64 mod = 13 * d;
        total += (i64)mus[i] * count_cong(L, R, mod, rem);
    }
    return total;
}

/* ---- binary search for max a in negative branch ---- */
static i64 max_abs_p_negative(i64 N) {
    i64 hi = isqrt64(N) + 2;
    i64 lo = 0;
    while (lo + 1 < hi) {
        i64 a = (lo + hi) / 2;
        if (a == 0) { lo = a; continue; }

        i64 qmin = (i64)(SQRT3 * a) + 1;
        i64 thr = 3 * a * a;
        while (qmin * qmin <= thr) qmin++;

        i64 qmax = (5 * a - 1) / 2;
        int ok;
        if (qmin > qmax) {
            ok = 0;
        } else {
            i64 z = -(qmin * qmin - 5 * a * qmin + 3 * a * a);
            ok = (z <= N);
        }
        if (ok) lo = a;
        else    hi = a;
    }
    return lo;
}

/* ---- main count C(N) ---- */
static i64 C_func(i64 N) {
    i64 fourN = 4 * N;
    i64 total = 0;

    int primes[32];
    i64 ds[256];
    int mus[256];

    /* Branch 1: p > 0, q > sqrt(3)*p, z = q^2 + 5pq + 3p^2 */
    i64 pmax = isqrt64(N / 3);
    for (i64 p = 1; p <= pmax; p++) {
        i64 qmin = (i64)(SQRT3 * p) + 1;
        i64 thr = 3 * p * p;
        while (qmin * qmin <= thr) qmin++;

        i64 disc = 13 * p * p + fourN;
        i64 qmax = (isqrt64(disc) - 5 * p) / 2;
        if (qmax < qmin) continue;

        i64 pp3 = 3 * p * p;
        while (qmax >= qmin && (qmax * qmax + 5 * p * qmax + pp3) > N) qmax--;
        if (qmax < qmin) continue;

        int np = distinct_prime_factors(p, primes);
        int nd = gen_divisors_mu(primes, np, ds, mus);

        i64 cnt = count_coprime_interval(qmin, qmax, ds, mus, nd);

        if (p % 13 != 0) {
            i64 bad = count_coprime_mod13(ds, mus, nd, qmin, qmax,
                                          (int)((4 * p) % 13));
            cnt -= bad;
        }
        total += cnt;
    }

    /* Branch 2: p = -a < 0, sqrt(3)*a < q < 2.5a, z = -(q^2 - 5aq + 3a^2) */
    i64 amax = max_abs_p_negative(N);
    i64 threshold = isqrt64(fourN / 13);

    for (i64 a = 1; a <= amax; a++) {
        i64 qmin = (i64)(SQRT3 * a) + 1;
        i64 thr = 3 * a * a;
        while (qmin * qmin <= thr) qmin++;

        i64 qmax = (5 * a - 1) / 2;
        if (qmax < qmin) continue;

        if (a > threshold) {
            i64 disc2 = 13 * a * a - fourN;
            i64 s = isqrt64(disc2);
            i64 lim = (5 * a - s) / 2;
            if (lim < qmax) qmax = lim;
        }

        while (qmax >= qmin &&
               (-(qmax * qmax - 5 * a * qmax + 3 * a * a)) > N) qmax--;
        if (qmax < qmin) continue;

        int np = distinct_prime_factors(a, primes);
        int nd = gen_divisors_mu(primes, np, ds, mus);

        i64 cnt = count_coprime_interval(qmin, qmax, ds, mus, nd);

        if (a % 13 != 0) {
            int rem13 = (int)(((-4 * a) % 13 + 13) % 13);
            i64 bad = count_coprime_mod13(ds, mus, nd, qmin, qmax, rem13);
            cnt -= bad;
        }
        total += cnt;
    }

    return total;
}

long long p769_native(void) {
    init_inv_mod13();

    i64 N = 100000000000000LL;  /* 10^14 */
    i64 pmax = isqrt64(N / 3);
    i64 amax = max_abs_p_negative(N);
    i64 maxp = pmax > amax ? pmax : amax;

    build_spf((int)maxp);

    i64 ans = C_func(N);

    free(spf_arr);
    return ans;
}
