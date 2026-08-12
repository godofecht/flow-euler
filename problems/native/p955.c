/* Project Euler 955: Finding Triangles
 *
 * Index of the 70th triangular term in the sequence.
 * Port of the Python reference solver.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef __int128 i128;
typedef long long i64;
typedef unsigned long long u64;
typedef unsigned __int128 u128;

/* gcd */
static i128 gcd128(i128 a, i128 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { i128 t = a % b; a = b; b = t; }
    return a;
}

/* (a * b) mod n for i128, using __int128 when safe, otherwise splitting.
 * Since n can be up to ~2^80, a*b can exceed __int128 range.
 * We use the fact that a, b < n < 2^127, so we split into 64-bit halves. */
static i128 mulmod128(i128 a, i128 b, i128 n) {
    /* If both fit in 63 bits, direct multiply fits in __int128 */
    if (a >= -((i128)1 << 62) && a < ((i128)1 << 62) &&
        b >= -((i128)1 << 62) && b < ((i128)1 << 62)) {
        return (i128)((u128)a * (u128)b % (u128)n);
    }
    /* Russian peasant / binary multiplication */
    a %= n;
    if (a < 0) a += n;
    b %= n;
    if (b < 0) b += n;
    i128 result = 0;
    while (b > 0) {
        if (b & 1) {
            result += a;
            if (result >= n) result -= n;
        }
        a <<= 1;
        if (a >= n) a -= n;
        b >>= 1;
    }
    return result;
}

static i128 powmod128(i128 base, i128 exp, i128 n) {
    i128 r = 1 % n, b = base % n;
    if (b < 0) b += n;
    while (exp > 0) {
        if (exp & 1) r = mulmod128(r, b, n);
        b = mulmod128(b, b, n);
        exp >>= 1;
    }
    return r;
}

static int is_prime128(i128 n) {
    if (n < 2) return 0;
    i64 small[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
    for (int i = 0; i < 10; i++) {
        if (n % small[i] == 0) return n == small[i];
    }
    i128 d = n - 1;
    i64 s = 0;
    while ((d & 1) == 0) { s++; d >>= 1; }
    i64 bases[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    for (int i = 0; i < 7; i++) {
        i128 a = (i128)bases[i] % n;
        if (a == 0) continue;
        i128 x = powmod128(a, d, n);
        if (x == 1 || x == n - 1) continue;
        int composite = 1;
        for (i64 j = 0; j < s - 1; j++) {
            x = mulmod128(x, x, n);
            if (x == n - 1) { composite = 0; break; }
        }
        if (composite) return 0;
    }
    return 1;
}

static i64 rng_state = 123456789;

static i64 rng_next(i64 modv) {
    rng_state = rng_state * 6364136223846793005LL + 1;
    i64 v = rng_state;
    if (v < 0) v = -v;
    return v % modv;
}

static i128 pollard_rho(i128 n) {
    if (n % 2 == 0) return 2;
    if (n % 3 == 0) return 3;
    while (1) {
        /* Use rng to get deterministic results */
        i128 x = (i128)(2 + rng_next(1000000007LL)) % n;
        i128 y = x;
        i128 c = (i128)(1 + rng_next(1000000007LL)) % n;
        i128 d = 1;
        while (d == 1) {
            x = (mulmod128(x, x, n) + c) % n;
            y = (mulmod128(y, y, n) + c) % n;
            y = (mulmod128(y, y, n) + c) % n;
            i128 diff = x - y;
            if (diff < 0) diff = -diff;
            d = gcd128(diff, n);
        }
        if (d != n) return d;
    }
}

typedef struct { i128 *p; i64 *e; i64 count; } Factors;

static void factor_rec(i128 n, Factors *f) {
    if (n == 1) return;
    if (is_prime128(n)) {
        for (i64 i = 0; i < f->count; i++) {
            if (f->p[i] == n) { f->e[i]++; return; }
        }
        f->p[f->count] = n;
        f->e[f->count] = 1;
        f->count++;
        return;
    }
    i128 d = pollard_rho(n);
    factor_rec(d, f);
    factor_rec(n / d, f);
}

static int cmp_i128(const void *a, const void *b) {
    i128 va = *(const i128 *)a;
    i128 vb = *(const i128 *)b;
    if (va < vb) return -1;
    if (va > vb) return 1;
    return 0;
}

static i128 *all_divisors(i128 n, i64 *out_count) {
    Factors f;
    f.p = malloc(64 * sizeof(i128));
    f.e = malloc(64 * sizeof(i64));
    f.count = 0;
    factor_rec(n, &f);

    i128 *ds = malloc(50000 * sizeof(i128));
    i64 nd = 1;
    ds[0] = 1;
    for (i64 i = 0; i < f.count; i++) {
        i64 old = nd;
        i128 pe = 1;
        for (i64 e = 0; e < f.e[i]; e++) {
            pe *= f.p[i];
            for (i64 j = 0; j < old; j++) {
                ds[nd++] = ds[j] * pe;
            }
        }
    }
    qsort(ds, nd, sizeof(i128), cmp_i128);
    *out_count = nd;
    free(f.p);
    free(f.e);
    return ds;
}

static void next_triangle_jump(i128 a_tri, i128 *out_k, i128 *out_a) {
    i128 M = 2 * a_tri;
    i64 nd;
    i128 *ds = all_divisors(M, &nd);
    for (i64 i = 0; i < nd; i++) {
        i128 x = ds[i];
        i128 y = M / x;
        if ((x - y) & 1) {
            i128 k = (x - y - 1) / 2;
            if (k > 0) {
                *out_k = k;
                *out_a = a_tri + k * (k + 1) / 2;
                free(ds);
                return;
            }
        }
    }
    free(ds);
    *out_k = -1;
    *out_a = -1;
}

long long p955_native(void) {
    i128 n = 0;
    i128 a = 3;
    for (int i = 0; i < 69; i++) {
        i128 step, na;
        next_triangle_jump(a, &step, &na);
        n += step;
        a = na;
    }
    return (long long)n;
}
