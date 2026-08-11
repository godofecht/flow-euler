/*
 * Project Euler 785: Symmetric Diophantine Equation
 *
 * S(N) = sum (x+y+z) over primitive solutions of
 *     15(x^2+y^2+z^2) = 34(xy+yz+zx)
 * with 1 <= x <= y <= z <= N and gcd(x,y,z)=1.
 *
 * Parameterisation into coprime pairs (a,b) with b < a, b < 3a/5.
 * The three coordinates are:
 *     A = 2ab + 3b^2
 *     B = 5a^2 - 2ab
 *     C = 3a^2 - 8ab + 5b^2
 * Contribution per valid pair is 8*(a^2 - ab + b^2).
 * Non-primitive exactly when all of A,B,C divisible by 19, which
 * reduces to A%19==0 and B%19==0.
 */

#include <stdint.h>
#include <stdlib.h>
#include <math.h>

typedef long long i64;
typedef unsigned long long u64;

/* Smallest prime factor sieve up to max_a. */
static int *spf_build(i64 max_a) {
    int *spf = (int *)malloc((size_t)(max_a + 1) * sizeof(int));
    if (!spf) return NULL;
    for (i64 i = 0; i <= max_a; i++) spf[i] = (int)i;
    for (i64 i = 2; i * i <= max_a; i++) {
        if (spf[i] == i) {
            for (i64 j = i * i; j <= max_a; j += i) {
                if (spf[j] == j) spf[j] = (int)i;
            }
        }
    }
    return spf;
}

static i64 isqrt_i64(i64 x) {
    if (x <= 0) return 0;
    i64 r = (i64)sqrt((double)x);
    /* adjust for floating error */
    while (r * r > x) r--;
    while ((r + 1) * (r + 1) <= x) r++;
    return r;
}

static i64 solve(i64 N) {
    /* a <= floor(sqrt(N/3)) since B = 5a^2 - 2ab >= 3a^2 (b < a). */
    i64 max_a = isqrt_i64(N / 3);
    if (max_a < 1) return 0;

    int *spf = spf_build(max_a);
    if (!spf) return 0;

    i64 N3 = 12 * N;
    i64 N5 = 5 * N;

    u64 total = 0;

    /* buffer for distinct prime factors of a */
    int pf[64];
    int pf_n;

    for (i64 a = 1; a <= max_a; a++) {
        i64 aa = a * a;

        /* distinct prime factors of a */
        pf_n = 0;
        {
            i64 x = a;
            while (x > 1) {
                int p = spf[x];
                pf[pf_n++] = p;
                while (x % p == 0) x /= p;
            }
        }

        /* bmin from B = 5a^2 - 2ab <= N  ==>  b >= (5a^2 - N)/(2a) */
        i64 bmin;
        if (5 * aa <= N) {
            bmin = 1;
        } else {
            i64 num = 5 * aa - N;
            bmin = (num + (2 * a - 1)) / (2 * a);
            if (bmin < 1) bmin = 1;
        }

        /* C > 0 requires b < 3a/5 */
        i64 bmax_pos = (3 * a - 1) / 5;
        if (bmax_pos < bmin) continue;

        /* A = 2ab + 3b^2 <= N  ->  b <= floor((-2a + sqrt(4a^2+12N))/6) */
        i64 disc1 = 4 * aa + N3;
        i64 bmax1 = (-2 * a + isqrt_i64(disc1)) / 6;

        /* C = 3a^2 - 8ab + 5b^2 <= N
         * 5b^2 - 8ab + (3a^2 - N) <= 0
         * roots: (4a +/- sqrt(a^2 + 5N))/5 */
        i64 disc3 = aa + N5;
        i64 sdisc3 = isqrt_i64(disc3);
        i64 bmax3 = (4 * a + sdisc3) / 5;
        i64 bmin3 = (4 * a - sdisc3 + 4) / 5; /* ceil */

        if (bmin3 > bmin) bmin = bmin3;

        i64 bmax = bmax_pos;
        if (bmax1 < bmax) bmax = bmax1;
        if (bmax3 < bmax) bmax = bmax3;
        if (bmax >= a) bmax = a - 1;

        if (bmax < bmin) continue;

        for (i64 b = bmin; b <= bmax; b++) {
            /* coprimality: b not divisible by any prime factor of a */
            int ok = 1;
            for (int k = 0; k < pf_n; k++) {
                if (b % pf[k] == 0) { ok = 0; break; }
            }
            if (!ok) continue;

            i64 ab = a * b;
            i64 bb = b * b;

            i64 A = 2 * ab + 3 * bb;
            i64 B = 5 * aa - 2 * ab;
            i64 C = 3 * aa - 8 * ab + 5 * bb;

            if (C <= 0) continue;
            if (A > N || B > N || C > N) continue;

            /* non-primitive iff all divisible by 19 */
            if (A % 19 == 0 && B % 19 == 0) continue;

            total += (u64)(8 * (aa - ab + bb));
        }
    }

    free(spf);
    return (i64)total;
}

long long p785_native(void) {
    return solve(1000000000LL);
}
