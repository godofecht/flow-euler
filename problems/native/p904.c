// Project Euler 904: Pythagorean Angle
// F(N, L) = sum_{n=1..N} f(cuberoot(n), L)
// N = 45000, L = 10^10
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define _DEG2RAD (M_PI / 180.0)
#define _T0 (sqrt(2.0) - 1.0)

static double g_of_t(double t) {
    double tt = t * t;
    double denom = 1.0 + tt;
    return 3.0 * t * (1.0 - tt) / (denom * denom);
}

static double root_left(double y) {
    double lo = 0.0, hi = _T0;
    for (int i = 0; i < 35; i++) {
        double mid = (lo + hi) * 0.5;
        if (g_of_t(mid) < y) lo = mid;
        else hi = mid;
    }
    return (lo + hi) * 0.5;
}

static double root_right(double y) {
    double lo = _T0, hi = 1.0;
    for (int i = 0; i < 35; i++) {
        double mid = (lo + hi) * 0.5;
        if (g_of_t(mid) > y) lo = mid;
        else hi = mid;
    }
    return (lo + hi) * 0.5;
}

typedef struct {
    long long p, q;
} Pair;

// Continued fraction candidates for real x in (0,1), p^2+q^2 <= L
static int cf_candidates_circle(double x, long long L, Pair *out, int max_out) {
    int count = 0;
    long long p0 = 0, q0 = 1, p1 = 1, q1 = 0;
    double frac = x;

    for (int iter = 0; iter < 80; iter++) {
        long long a = (long long)frac;
        // Check for potential overflow before computing p2, q2
        // p2 = a * p1 + p0, q2 = a * q1 + q0
        // If a is huge, p2/q2 will be huge and certainly > sqrt(L)
        long long sqrtL = (long long)sqrt((double)L) + 1;
        long long p2, q2;
        if (a > 2 * sqrtL || (p1 > 0 && a > sqrtL / p1 + 1)) {
            // p2 will certainly exceed sqrt(L), so p2^2+q2^2 > L
            // Treat as if p2*p2+q2*q2 > L
            p2 = 0; q2 = 0; // dummy, won't be used
            // Go to semiconvergent branch
            long long A = p1 * p1 + q1 * q1;
            long long B = 2 * (p0 * p1 + q0 * q1);
            long long C = p0 * p0 + q0 * q0 - L;

            long long kmax = 0;
            if (A > 0) {
                __int128 disc = (__int128)B * B - (__int128)4 * A * C;
                if (disc > 0) {
                    long long s_val = (long long)sqrt((double)disc);
                    while ((__int128)s_val * s_val > disc) s_val--;
                    while ((__int128)(s_val + 1) * (s_val + 1) <= disc) s_val++;
                    kmax = (-B + s_val) / (2 * A);
                }
            }
            if (kmax > a - 1) kmax = a - 1;
            if (kmax < 0) kmax = 0;

            long long lo = kmax - 3; if (lo < 1) lo = 1;
            long long hi = kmax + 3; if (hi > a - 1) hi = a - 1;
            if (hi < lo) hi = lo; // ensure valid range

            for (long long k = lo; k <= hi; k++) {
                long long ps = k * p1 + p0;
                long long qs = k * q1 + q0;
                if (ps > 0 && qs > 0 && (__int128)ps * ps + (__int128)qs * qs <= L) {
                    if (count < max_out) {
                        out[count].p = ps;
                        out[count].q = qs;
                        count++;
                    }
                }
            }
            break;
        }

        p2 = a * p1 + p0;
        q2 = a * q1 + q0;

        if ((__int128)p2 * p2 + (__int128)q2 * q2 > L) {
            // Find largest k <= a-1 with (p0+k*p1)^2 + (q0+k*q1)^2 <= L
            long long A = p1 * p1 + q1 * q1;
            long long B = 2 * (p0 * p1 + q0 * q1);
            long long C = p0 * p0 + q0 * q0 - L;

            long long kmax = 0;
            if (A > 0) {
                __int128 disc = (__int128)B * B - (__int128)4 * A * C;
                if (disc > 0) {
                    long long s_val = (long long)sqrt((double)disc);
                    while ((__int128)s_val * s_val > disc) s_val--;
                    while ((__int128)(s_val + 1) * (s_val + 1) <= disc) s_val++;
                    kmax = (-B + s_val) / (2 * A);
                }
            }
            if (kmax > a - 1) kmax = a - 1;
            if (kmax < 0) kmax = 0;

            long long lo = kmax - 3; if (lo < 1) lo = 1;
            long long hi = kmax + 3; if (hi > a - 1) hi = a - 1;

            for (long long k = lo; k <= hi; k++) {
                long long ps = k * p1 + p0;
                long long qs = k * q1 + q0;
                if (ps > 0 && qs > 0 && (__int128)ps * ps + (__int128)qs * qs <= L) {
                    if (count < max_out) {
                        out[count].p = ps;
                        out[count].q = qs;
                        count++;
                    }
                }
            }
            break;
        }

        if (p2 > 0 && q2 > 0) {
            if (count < max_out) {
                out[count].p = p2;
                out[count].q = q2;
                count++;
            }
        }

        if (frac == (double)a) break;
        frac = 1.0 / (frac - a);
        p0 = p1; q0 = q1; p1 = p2; q1 = q2;
    }

    // Also include last convergent if admissible
    if (p1 > 0 && q1 > 0 && (__int128)p1 * p1 + (__int128)q1 * q1 <= L) {
        int found = 0;
        for (int i = 0; i < count; i++) {
            if (out[i].p == p1 && out[i].q == q1) { found = 1; break; }
        }
        if (!found && count < max_out) {
            out[count].p = p1;
            out[count].q = q1;
            count++;
        }
    }

    return count;
}

static void triangle_from_mn(long long m, long long n, long long *a, long long *b, long long *c) {
    long long aa = m * m - n * n;
    long long bb = 2LL * m * n;
    long long cc = m * m + n * n;
    if ((m + n) % 2 == 0) {
        *a = aa / 2; *b = bb / 2; *c = cc / 2;
    } else {
        *a = aa; *b = bb; *c = cc;
    }
}

static double tan_theta_from_legs(long long a, long long b) {
    double aa = (double)a * a;
    double bb = (double)b * b;
    return (3.0 * a * b) / (2.0 * (aa + bb));
}

static long long f_val(double alpha_deg, long long L) {
    double alpha_rad = alpha_deg * _DEG2RAD;
    double y = tan(alpha_rad);

    double r1 = root_left(y);
    double r2 = root_right(y);

    double best_diff = 1e300;
    double best_area_key = -1;
    long long best_perim = 0;

    double roots[2] = {r1, r2};
    Pair cands[200];

    for (int ri = 0; ri < 2; ri++) {
        int ncands = cf_candidates_circle(roots[ri], L, cands, 200);
        for (int ci = 0; ci < ncands; ci++) {
            long long n_val = cands[ci].p;  // p is the smaller (n)
            long long m_val = cands[ci].q;  // q is the larger (m)
            if (!(0 < n_val && n_val < m_val)) continue;

            long long a, b, c;
            triangle_from_mn(m_val, n_val, &a, &b, &c);
            if (a <= 0 || b <= 0) continue;
            if (c > L) continue;

            long long k = L / c;
            if (k <= 0) continue;

            double theta = atan(tan_theta_from_legs(a, b));
            double diff = fabs(theta - alpha_rad);

            double area_key = (double)(k * k) * (double)a * b;
            long long perim = k * (a + b + c);

            if (diff + 1e-16 < best_diff) {
                best_diff = diff;
                best_area_key = area_key;
                best_perim = perim;
            } else if (fabs(diff - best_diff) <= 1e-16) {
                if (area_key > best_area_key) {
                    best_area_key = area_key;
                    best_perim = perim;
                }
            }
        }
    }

    return best_perim;
}

long long p904_native(void) {
    long long N = 45000;
    long long L = 10000000000LL;  // 10^10
    long long total = 0;
    for (long long n = 1; n <= N; n++) {
        double alpha = cbrt((double)n);
        total += f_val(alpha, L);
    }
    return total;
}
