// Project Euler 891: Ambiguous Clock
// Count ambiguous moments in a 12-hour cycle with three identical hands.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long i64;

#define L 43200
#define SHIFT (1LL << 20)

static i64 gcd(i64 a, i64 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { i64 t = a % b; a = b; b = t; }
    return a;
}

/* ceil(n/d) for d>0, matching Python's -((-n)//d) with floor division */
static i64 ceil_div(i64 n, i64 d) {
    /* d > 0 */
    if (n >= 0) return (n + d - 1) / d;
    return -((-n) / d);
}

/* floor(n/d) for d>0, matching Python's n // d */
static i64 floor_div(i64 n, i64 d) {
    /* d > 0 */
    if (n >= 0) return n / d;
    return -((-n + d - 1) / d);
}

static void k2_range(i64 u, i64 v, i64 D, i64 *lo, i64 *hi) {
    if (v > 0) {
        *lo = ceil_div(u - D + 1, v);
        *hi = floor_div(u, v);
    } else {
        i64 vp = -v;
        *lo = ceil_div(-u, vp);
        *hi = floor_div(D - 1 - u, vp);
    }
}

static int cmp_i64(const void *a, const void *b) {
    i64 x = *(const i64 *)a, y = *(const i64 *)b;
    if (x < y) return -1;
    if (x > y) return 1;
    return 0;
}

long long p891_native(void) {
    static const i64 V[3] = {1, 12, 720};
    i64 a = V[0] - V[1]; /* -11 */
    i64 c = V[0] - V[2]; /* -719 */

    int perms[5][3] = {
        {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}
    };

    /* Collect all keys, then sort and count unique */
    i64 cap = 5000000;
    i64 *keys = (i64 *)malloc(cap * sizeof(i64));
    i64 nkeys = 0;

    for (int pi = 0; pi < 5; pi++) {
        int s0 = perms[pi][0], s1 = perms[pi][1], s2 = perms[pi][2];
        i64 b = -(V[s0] - V[s1]);
        i64 d = -(V[s0] - V[s2]);
        i64 det = a * d - c * b;
        i64 sign = 1;
        if (det < 0) { sign = -1; det = -det; }
        i64 D = det;
        i64 K1 = (a < 0 ? -a : a) + (b < 0 ? -b : b);
        i64 K2 = (c < 0 ? -c : c) + (d < 0 ? -d : d);

        i64 g0 = gcd(L, D);
        i64 l1 = L / g0;
        i64 d1 = D / g0;

        for (i64 k1 = -K1; k1 <= K1; k1++) {
            i64 u1 = sign * k1 * d;
            i64 v1 = sign * b;
            i64 u2 = sign * 719 * k1;
            i64 v2 = sign * 11;

            i64 lo1, hi1, lo2, hi2;
            k2_range(u1, v1, D, &lo1, &hi1);
            k2_range(u2, v2, D, &lo2, &hi2);

            i64 lo = -K2 > lo1 ? -K2 : lo1;
            if (lo2 > lo) lo = lo2;
            i64 hi = K2 < hi1 ? K2 : hi1;
            if (hi2 < hi) hi = hi2;

            if (lo > hi) continue;

            for (i64 k2 = lo; k2 <= hi; k2++) {
                i64 num_t = u1 - v1 * k2;
                i64 num_tp = u2 - v2 * k2;
                if (num_t == num_tp) continue;

                i64 g = gcd(num_t, d1);
                i64 n = l1 * (num_t / g);
                i64 den = d1 / g;
                keys[nkeys++] = n * SHIFT + den;

                g = gcd(num_tp, d1);
                n = l1 * (num_tp / g);
                den = d1 / g;
                keys[nkeys++] = n * SHIFT + den;
            }
        }
    }

    qsort(keys, nkeys, sizeof(i64), cmp_i64);

    i64 unique = 0;
    for (i64 i = 0; i < nkeys; i++) {
        if (i == 0 || keys[i] != keys[i - 1]) unique++;
    }

    free(keys);
    return unique;
}
