// Project Euler 975
// Port of the Python reference solver.
// Winding walk on turning points of H(a,b,x) over prime pairs.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { int num, den; } Point;

static int gcd_int(int a, int b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int t = a % b; a = b; b = t; }
    return a;
}

static Point make_normalized(int num, int den) {
    int g = gcd_int(num, den);
    Point p = { num / g, den / g };
    return p;
}

static double height_val(int a, int b, int num, int den) {
    if (num == 0) return 0.0;
    if (num == den) return 1.0;
    double x = (double)num / (double)den;
    double z = 0.5 - (b * cos((double)a * M_PI * x) + a * cos((double)b * M_PI * x)) / (2.0 * (a + b));
    if (z < 0.0 && z > -1e-14) return 0.0;
    if (z > 1.0 && z < 1.0 + 1e-14) return 1.0;
    return z;
}

static int deriv_sign(int a, int b, int ln, int ld, int rn, int rd) {
    double x = (double)(ln * rd + rn * ld) / (2.0 * (double)ld * (double)rd);
    double value = sin((double)(a + b) * M_PI * x * 0.5) * cos((double)abs(a - b) * M_PI * x * 0.5);
    return (value > 0.0) ? 1 : -1;
}

static int cmp_points(const void *pa, const void *pb) {
    const Point *a = (const Point *)pa;
    const Point *b = (const Point *)pb;
    double va = (double)a->num / (double)a->den;
    double vb = (double)b->num / (double)b->den;
    if (va < vb) return -1;
    if (va > vb) return 1;
    return 0;
}

#define MAX_CAND 4096

static Point g_cand[MAX_CAND];
static int g_signs[MAX_CAND];
static Point g_kept[MAX_CAND];

static int turning_values(int a, int b, double *tv) {
    int s = a + b;
    int delta = abs(a - b);

    int nc = 0;
    for (int k = 0; k <= s / 2; k++)
        g_cand[nc++] = make_normalized(2 * k, s);
    for (int k = 0; k < delta / 2; k++)
        g_cand[nc++] = make_normalized(2 * k + 1, delta);

    qsort(g_cand, nc, sizeof(Point), cmp_points);

    int unique = 1;
    for (int i = 1; i < nc; i++) {
        if (g_cand[i].num != g_cand[unique - 1].num ||
            g_cand[i].den != g_cand[unique - 1].den) {
            g_cand[unique++] = g_cand[i];
        }
    }
    nc = unique;

    for (int i = 0; i < nc - 1; i++) {
        g_signs[i] = deriv_sign(a, b,
            g_cand[i].num, g_cand[i].den,
            g_cand[i + 1].num, g_cand[i + 1].den);
    }

    int nk = 0;
    g_kept[nk++] = g_cand[0];
    for (int i = 1; i < nc - 1; i++) {
        if (g_signs[i - 1] != g_signs[i])
            g_kept[nk++] = g_cand[i];
    }
    g_kept[nk++] = g_cand[nc - 1];

    for (int i = 0; i < nk; i++)
        tv[i] = height_val(a, b, g_kept[i].num, g_kept[i].den);

    return nk;
}

static double F_func(int a, int b, int c, int d) {
    double za[MAX_CAND], zb[MAX_CAND];
    int na = turning_values(a, b, za);
    int nb = turning_values(c, d, zb);

    int i = 0, j = 0;
    double current = 0.0, total = 0.0;
    int upward = 1;
    double eps = 1e-12;
    long long max_steps = 4LL * (na + nb) * (na + nb);

    for (long long step = 0; step < max_steps; step++) {
        if (i < 0 || i >= na - 1 || j < 0 || j >= nb - 1) {
            if (fabs(current - 1.0) < 1e-9) return total;
            return total;
        }
        double a0 = za[i], a1 = za[i + 1];
        double b0 = zb[j], b1 = zb[j + 1];
        double lower = fmax(fmin(a0, a1), fmin(b0, b1));
        double upper = fmin(fmax(a0, a1), fmax(b0, b1));
        double nxt = upward ? upper : lower;
        total += fabs(nxt - current);
        int advanced = 0;
        if (fabs(nxt - a0) <= eps) { i -= 1; advanced = 1; }
        else if (fabs(nxt - a1) <= eps) { i += 1; advanced = 1; }
        if (fabs(nxt - b0) <= eps) { j -= 1; advanced = 1; }
        else if (fabs(nxt - b1) <= eps) { j += 1; advanced = 1; }
        if (!advanced) return total;
        current = nxt;
        upward = !upward;
    }
    return total;
}

static double G_func(int m, int n) {
    char *is_prime = (char *)malloc(n + 1);
    memset(is_prime, 1, n + 1);
    is_prime[0] = is_prime[1] = 0;
    for (int p = 2; (long long)p * p <= n; p++) {
        if (is_prime[p]) {
            for (int j = p * p; j <= n; j += p) is_prime[j] = 0;
        }
    }

    int ps[256];
    int np = 0;
    for (int p = m; p <= n; p++) {
        if (is_prime[p]) ps[np++] = p;
    }
    free(is_prime);

    double total = 0.0;
    for (int i = 0; i < np; i++) {
        for (int j = i + 1; j < np; j++) {
            int p = ps[i], q = ps[j];
            total += F_func(p, q, p, 2 * q - p);
        }
    }
    return total;
}

double p975_native(void) {
    return G_func(500, 1000);
}
