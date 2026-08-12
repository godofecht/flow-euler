/* Project Euler 950: Pirate Treasure
 * Compute sum_{k=1..6} T(10^16, 10^k+1, 1/sqrt(10^k+1)) mod 10^9.
 * All arithmetic is exact integer (no floating point).
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

typedef int64_t i64;
typedef __int128 i128;

#define MOD_9 1000000000LL

/* Integer square root of a __int128 value */
static i128 isqrt128(i128 x) {
    if (x <= 0) return 0;
    if (x < (i128)1 << 63) {
        /* fits in i64, use hardware-friendly path */
        i64 r = (i64)x;
        i64 s = 0;
        i64 b = 1LL << 62;
        while (b > r) b >>= 2;
        while (b > 0) {
            if (r >= s + b) { r -= s + b; s = (s >> 1) + b; }
            else { s >>= 1; }
            b >>= 2;
        }
        return s;
    }
    /* Newton's method for large values */
    i128 r = (i128)1 << 63; /* initial guess: 2^63 */
    /* Find a reasonable starting point */
    int bits = 0;
    i128 tmp = x;
    while (tmp) { bits++; tmp >>= 1; }
    r = (i128)1 << (bits / 2 + 1);
    while (1) {
        i128 next = (r + x / r) / 2;
        if (next >= r) break;
        r = next;
    }
    /* Correct: r might be off by 1 */
    while (r * r > x) r--;
    while ((r + 1) * (r + 1) <= x) r++;
    return r;
}

/* floor(d / sqrt(D)) for d >= 0, D >= 1 */
static i128 floor_div_sqrt(i128 d, i128 D) {
    if (d <= 0) return 0;
    i128 dd = d * d;
    i128 t = isqrt128(dd / D);
    while ((t + 1) * (t + 1) * D <= dd) t++;
    while (t * t * D > dd) t--;
    return t;
}

/* ceil(d / sqrt(D)) for d >= 0, D non-square */
static i128 ceil_div_sqrt(i128 d, i128 D) {
    if (d <= 0) return 0;
    return floor_div_sqrt(d, D) + 1;
}

static i128 initial_prefix_sum(i128 N, i128 C) {
    if (N <= 0) return 0;
    i128 limit = 2 * C + 2;
    i128 M = N < limit ? N : limit;
    if (M > 2 * C) M = 2 * C;
    if (M <= 0) return 0;

    i128 m = M / 2;
    i128 s = 2 * (m * C - (m * (m - 1)) / 2);
    if (M % 2 == 1) {
        s += C - m;
    }
    return s;
}

static i128 next_reset(i128 L, i128 C, i128 D) {
    if (C == 0) {
        return 2 * L;
    }

    i128 t = 1;
    while (t <= C) {
        i128 y = C / t;
        i128 x = 2 * L - 2 * y;
        i128 d = x - L;
        if (d > 0) {
            i128 s = ceil_div_sqrt(d, D);
            if (C / s == y) {
                return x;
            }
        }
        t = C / y + 1;
    }

    /* y = 0 case (k > C) */
    return 2 * L;
}

static i128 T_func(i128 N, i128 C, i128 D) {
    if (N <= 0) return 0;

    i128 start_reset = 2 * C + 2;
    if (N <= start_reset) {
        return initial_prefix_sum(N, C);
    }

    i128 total = initial_prefix_sum(start_reset, C);
    i128 L = start_reset;
    i128 cL = 0;

    while (L < N) {
        i128 x = next_reset(L, C, D);
        if (x > N) {
            i128 d = N - L + 1;
            total += (d - 1) * cL + (d - 1) * d / 2;
            break;
        }

        i128 d = x - L;
        if (d > 1) {
            total += (d - 1) * cL + (d - 1) * d / 2;
        }

        i128 required_votes = (x + 1) / 2;
        i128 free_votes = d;
        i128 need_bribes = required_votes - free_votes;
        if (need_bribes < 0) need_bribes = 0;

        i128 s = ceil_div_sqrt(d, D);
        i128 cost = need_bribes * s;
        cL = C - cost;

        total += cL;
        L = x;
    }

    return total;
}

long long p950_native(void) {
    i128 N = (i128)1 << 0;
    /* N = 10^16 */
    i128 Nval = 1;
    for (int i = 0; i < 16; i++) Nval *= 10;

    i128 acc = 0;
    for (int k = 1; k <= 6; k++) {
        i128 C = 1;
        for (int i = 0; i < k; i++) C *= 10;
        C += 1;
        acc += T_func(Nval, C, C);
    }

    i128 result = acc % MOD_9;
    if (result < 0) result += MOD_9;
    return (long long)result;
}
