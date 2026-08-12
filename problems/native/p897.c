#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 128

/* Tridiagonal solve (Thomas algorithm) */
static void tridiag_solve(const double *lower, const double *diag,
                          const double *upper, const double *rhs,
                          double *x, int n) {
    double c[MAXN], d[MAXN], b[MAXN];
    memcpy(c, upper, (n - 1) * sizeof(double));
    memcpy(d, diag, n * sizeof(double));
    memcpy(b, rhs, n * sizeof(double));
    for (int i = 0; i < n - 1; i++) {
        double w = lower[i] / d[i];
        d[i + 1] -= w * c[i];
        b[i + 1] -= w * b[i];
    }
    x[n - 1] = b[n - 1] / d[n - 1];
    for (int i = n - 2; i >= 0; i--)
        x[i] = (b[i] - c[i] * x[i + 1]) / d[i];
}

static double max_residual(const double *xs, int n) {
    double mx = 0.0;
    for (int k = 1; k < n - 1; k++) {
        double a = xs[k - 1], x = xs[k], b = xs[k + 1];
        double a4 = a * a * a * a;
        double b4 = b * b * b * b;
        double x3 = x * x * x;
        double g = a4 - b4 + 4.0 * x3 * (b - a);
        double ag = fabs(g);
        if (ag > mx) mx = ag;
    }
    return mx;
}

static void newton_solve(double *xs, int n) {
    int m = n - 2;
    static double g[MAXN], diag[MAXN], lower[MAXN], upper[MAXN], rhs[MAXN], delta[MAXN];

    for (int iter = 0; iter < 200; iter++) {
        memset(g, 0, m * sizeof(double));
        memset(diag, 0, m * sizeof(double));
        memset(lower, 0, (m - 1) * sizeof(double));
        memset(upper, 0, (m - 1) * sizeof(double));

        for (int k = 1; k < n - 1; k++) {
            double a = xs[k - 1], x = xs[k], b = xs[k + 1];
            int idx = k - 1;
            double a4 = a * a * a * a;
            double b4 = b * b * b * b;
            double x3 = x * x * x;
            double x2 = x * x;
            double a3 = a * a * a;
            double b3 = b * b * b;

            g[idx] = a4 - b4 + 4.0 * x3 * (b - a);
            diag[idx] = 12.0 * x2 * (b - a);
            if (idx - 1 >= 0)
                lower[idx - 1] = 4.0 * a3 - 4.0 * x3;
            if (idx + 1 < m)
                upper[idx] = 4.0 * x3 - 4.0 * b3;
        }

        double maxg = 0.0;
        for (int i = 0; i < m; i++)
            if (fabs(g[i]) > maxg) maxg = fabs(g[i]);
        if (maxg < 1e-15) break;

        for (int i = 0; i < m; i++) rhs[i] = -g[i];
        tridiag_solve(lower, diag, upper, rhs, delta, m);

        double alpha = 1.0;
        int found = 0;
        static double trial[MAXN];
        while (alpha > 1e-14) {
            memcpy(trial, xs, n * sizeof(double));
            for (int i = 0; i < m; i++)
                trial[i + 1] = xs[i + 1] + alpha * delta[i];

            int ok = 1;
            for (int i = 0; i < n - 1; i++) {
                if (!(trial[i] < trial[i + 1])) { ok = 0; break; }
            }
            if (ok && trial[1] > -1.0 && trial[1] < 1.0 &&
                trial[n - 2] > -1.0 && trial[n - 2] < 1.0) {
                if (max_residual(trial, n) < maxg) {
                    memcpy(xs, trial, n * sizeof(double));
                    found = 1;
                    break;
                }
            }
            alpha *= 0.5;
        }
        if (!found) {
            for (int i = 0; i < m; i++)
                xs[i + 1] += 1e-16 * delta[i];
        }
    }
}

static double kahan_sum(const double *values, int n) {
    double s = 0.0, c = 0.0;
    for (int i = 0; i < n; i++) {
        double y = values[i] - c;
        double t = s + y;
        c = (t - s) - y;
        s = t;
    }
    return s;
}

static int cmp_double(const void *a, const void *b) {
    double va = *(const double *)a;
    double vb = *(const double *)b;
    if (va < vb) return -1;
    if (va > vb) return 1;
    return 0;
}

static double compute_area(const double *xs, int n) {
    static double terms[MAXN];
    for (int i = 0; i < n - 1; i++) {
        double a = xs[i], b = xs[i + 1];
        double a4 = a * a * a * a;
        double b4 = b * b * b * b;
        terms[i] = (b - a) * (a4 + b4) * 0.5;
    }
    double bottom = kahan_sum(terms, n - 1);
    return 2.0 - bottom;
}

double p897_native(void) {
    int n = 101;
    int m = (n - 1) / 2;
    double best_area = -1e18;

    for (int extra_left = 0; extra_left <= 1; extra_left++) {
        int left_internal = m;
        int right_internal = m - 1;
        if (!extra_left) {
            int tmp = left_internal;
            left_internal = right_internal;
            right_internal = tmp;
        }

        double xs[MAXN];
        int cnt = 0;
        xs[cnt++] = -1.0;

        for (int j = 1; j <= left_internal; j++) {
            double z = 1.0 - (double)j / (left_internal + 1);
            xs[cnt++] = -pow(z, 3.0 / 5.0);
        }
        for (int j = 1; j <= right_internal; j++) {
            double z = (double)j / (right_internal + 1);
            xs[cnt++] = pow(z, 3.0 / 5.0);
        }
        xs[cnt++] = 1.0;

        qsort(xs, cnt, sizeof(double), cmp_double);
        newton_solve(xs, n);

        double area = compute_area(xs, n);
        if (area > best_area) best_area = area;
    }

    return best_area;
}
