#include <stdint.h>
#include <math.h>

typedef long long i64;
typedef __int128 i128;

/* Project Euler 919 - Fortunate Triangles
   |cos(A)| = 1/4 for at least one angle.
   S(P) = sum of perimeters of all fortunate triangles with perimeter <= P.
   P = 10^7.
*/

static i64 gcd64(i64 a, i64 b) {
    while (b) { i64 t = a % b; a = b; b = t; }
    return a;
}

static int passes_d_filter(i64 m, i64 d) {
    if (d == 1) return (m % 3 != 0) && (m % 5 != 0);
    if (d == 3) return (m % 5 != 0);
    if (d == 5) return (m % 3 != 0);
    return 1; /* d == 15 */
}

static i128 total;

static void process(i64 raw_a, i64 raw_b, i64 raw_c, i64 P) {
    i64 g = gcd64(raw_a, raw_b);
    g = gcd64(g, raw_c);
    i64 a = raw_a / g;
    i64 b = raw_b / g;
    i64 c = raw_c / g;
    /* Avoid double-counting: skip if a > b */
    if (a > b) return;
    i64 per = a + b + c;
    if (per > P) return;
    i64 n = P / per;
    total += (i128)per * (i128)n * (i128)(n + 1) / 2;
}

long long p919_native(void) {
    i64 P = 10000000;
    total = 0;
    i64 q_max = 2 * (i64)sqrt((double)P) + 10;
    i64 ds[] = {1, 3, 5, 15};

    /* Family A: cos = +1/4
       A = 8*p*q, B = 15*q^2 - p^2 + 2*p*q, C = 15*q^2 + p^2
       Constraint: p < 5q, perimeter raw = 10*q*(p+3q)
    */
    for (i64 q = 1; q <= q_max; q++) {
        i64 qq = q * q;
        i64 max_p = 5 * q - 1;
        for (int di = 0; di < 4; di++) {
            i64 d = ds[di];
            i64 p_lim = (P * 8 * d) / (10 * q) - 3 * q;
            if (p_lim > max_p) p_lim = max_p;
            if (p_lim < d) continue;
            for (i64 p = d; p <= p_lim; p += d) {
                i64 m = p / d;
                if (!passes_d_filter(m, d)) continue;
                if (gcd64(p, q) != 1) continue;
                i64 a = 8 * p * q;
                i64 b = 15 * qq - p * p + 2 * p * q;
                i64 c = 15 * qq + p * p;
                if (b <= 0) continue;
                process(a, b, c, P);
            }
        }
    }

    /* Family B: cos = -1/4
       A = 8*p*q, B = 15*q^2 - p^2 - 2*p*q, C = 15*q^2 + p^2
       Constraint: p < 3q, perimeter raw = 6*q*(p+5q)
    */
    for (i64 q = 1; q <= q_max; q++) {
        i64 qq = q * q;
        i64 max_p = 3 * q - 1;
        for (int di = 0; di < 4; di++) {
            i64 d = ds[di];
            i64 p_lim = (P * 8 * d) / (6 * q) - 5 * q;
            if (p_lim > max_p) p_lim = max_p;
            if (p_lim < d) continue;
            for (i64 p = d; p <= p_lim; p += d) {
                i64 m = p / d;
                if (!passes_d_filter(m, d)) continue;
                if (gcd64(p, q) != 1) continue;
                i64 a = 8 * p * q;
                i64 b = 15 * qq - p * p - 2 * p * q;
                i64 c = 15 * qq + p * p;
                if (b <= 0) continue;
                process(a, b, c, P);
            }
        }
    }

    return (i64)total;
}
