
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double arith_sum(long a, long b) {
    if (a > b) return 0.0;
    double n = (double)(b - a + 1);
    return ((double)(a + b)) * n * 0.5;
}

static double expected_min_consecutive(long low, long high) {
    long L = high - low + 1;
    double denom = (double)(L * L);
    double total = 0.0;
    for (long i = 0; i < L; i++) {
        double v = (double)(low + i);
        double w = (double)(2 * L - 2 * i - 1) / denom;
        total += v * w;
    }
    return total;
}

static void prefix_fill(const double *arr, double *pref, long M) {
    double s = 0.0;
    pref[0] = 0.0;
    for (long i = 1; i <= M; i++) {
        s += arr[i];
        pref[i] = s;
    }
}

static int compute_C1_iterative(long n, long m, double *C1, double *B0_out) {
    long L = m - n + 1;
    long A_low = n + 5, A_high = m + 5, M = A_high;
    double EminA = expected_min_consecutive(A_low, A_high);
    double k = (double)L / (double)(L - 1);
    double invL = 1.0 / (double)L;
    double meanA = 0.5 * (double)(A_low + A_high);
    double *C0 = calloc(M + 1, sizeof(double));
    double *newC0 = calloc(M + 1, sizeof(double));
    double *newC1 = calloc(M + 1, sizeof(double));
    double *pref = calloc(M + 1, sizeof(double));
    for (long x = 1; x <= M; x++) {
        double mx = meanA < (double)x ? meanA : (double)x;
        C1[x] = meanA + mx;
        C0[x] = meanA + (double)x;
    }
    double w = 1.15;
    int ok = 0;
    for (int it = 0; it < 2000; it++) {
        prefix_fill(C1, pref, M);
        double ss = 0.0;
        for (long d = 1; d < L; d++) ss += (double)(L - d) * 2.0 * C1[d];
        double B0 = k * (EminA + ss / (double)(L * L));
        for (long x = 1; x <= M; x++) {
            long lt_low = A_low, lt_high = A_high;
            if (lt_high > x - 1) lt_high = x - 1;
            double sum_lt_a = arith_sum(lt_low, lt_high);
            double sum_lt_C = 0.0;
            if (lt_low <= lt_high) {
                long y_lo = x - lt_high, y_hi = x - lt_low;
                sum_lt_C = pref[y_hi] - (y_lo > 1 ? pref[y_lo - 1] : 0.0);
            }
            long gt_low = A_low; if (gt_low < x + 1) gt_low = x + 1;
            long gt_high = A_high;
            long count_gt = gt_high - gt_low + 1; if (count_gt < 0) count_gt = 0;
            double base_gt = (double)count_gt * (double)x;
            double sum_gt_C = 0.0;
            if (gt_low <= gt_high) {
                long y_lo = gt_low - x, y_hi = gt_high - x;
                sum_gt_C = pref[y_hi] - (y_lo > 1 ? pref[y_lo - 1] : 0.0);
            }
            double total = sum_lt_a + sum_lt_C + base_gt + sum_gt_C;
            if (A_low <= x && x <= A_high) total += (double)x + B0;
            newC0[x] = total * invL;
        }
        for (long x = 1; x <= M; x++) C0[x] += w * (newC0[x] - C0[x]);
        prefix_fill(C0, pref, M);
        ss = 0.0;
        for (long d = 1; d < L; d++) ss += (double)(L - d) * C0[d];
        double B1 = k * (EminA + ss / (double)(L * L));
        double err = 0.0;
        int bad = 0;
        for (long x = 1; x <= M; x++) {
            long lt_low = A_low, lt_high = A_high;
            if (lt_high > x - 1) lt_high = x - 1;
            double sum_lt_a = arith_sum(lt_low, lt_high);
            long gt_low = A_low; if (gt_low < x + 1) gt_low = x + 1;
            long gt_high = A_high;
            long count_gt = gt_high - gt_low + 1; if (count_gt < 0) count_gt = 0;
            double base_gt = (double)count_gt * (double)x;
            double sum_gt_C = 0.0;
            if (gt_low <= gt_high) {
                long y_lo = gt_low - x, y_hi = gt_high - x;
                sum_gt_C = pref[y_hi] - (y_lo > 1 ? pref[y_lo - 1] : 0.0);
            }
            double total = sum_lt_a + base_gt + sum_gt_C;
            if (A_low <= x && x <= A_high) total += (double)x + B1;
            newC1[x] = total * invL;
            double v = C1[x] + w * (newC1[x] - C1[x]);
            if (v != v || fabs(v) > 1e300) { bad = 1; break; }
            double e = fabs(v - C1[x]);
            if (e > err) err = e;
            C1[x] = v;
        }
        if (bad) break;
        if (err < 1e-12) {
            ss = 0.0;
            for (long d = 1; d < L; d++) ss += (double)(L - d) * 2.0 * C1[d];
            *B0_out = k * (EminA + ss / (double)(L * L));
            ok = 1;
            break;
        }
    }
    free(C0); free(newC0); free(newC1); free(pref);
    return ok;
}

static void gauss_solve(double *A, double *b, long n, double *xout) {
    for (long i = 0; i < n; i++) A[i * (n + 1) + n] = b[i];
    for (long col = 0; col < n; col++) {
        long piv = col;
        double best = fabs(A[col * (n + 1) + col]);
        for (long r = col + 1; r < n; r++) {
            double v = fabs(A[r * (n + 1) + col]);
            if (v > best) { best = v; piv = r; }
        }
        if (piv != col) {
            for (long j = col; j <= n; j++) {
                double tmp = A[col * (n + 1) + j];
                A[col * (n + 1) + j] = A[piv * (n + 1) + j];
                A[piv * (n + 1) + j] = tmp;
            }
        }
        double inv = 1.0 / A[col * (n + 1) + col];
        for (long j = col; j <= n; j++) A[col * (n + 1) + j] *= inv;
        for (long r = 0; r < n; r++) if (r != col) {
            double factor = A[r * (n + 1) + col];
            if (factor == 0.0) continue;
            for (long j = col; j <= n; j++) A[r * (n + 1) + j] -= factor * A[col * (n + 1) + j];
        }
    }
    for (long i = 0; i < n; i++) xout[i] = A[i * (n + 1) + n];
}

static void compute_C1_direct(long n, long m, double *C1, double *B0_out) {
    long L = m - n + 1;
    long A_low = n + 5, A_high = m + 5, M = A_high;
    double EminA = expected_min_consecutive(A_low, A_high);
    double k = (double)L / (double)(L - 1);
    double invL = 1.0 / (double)L;
    double *q0 = calloc(M + 1, sizeof(double));
    double *p0 = calloc((M + 1) * (M + 1), sizeof(double));
    double *pref_q = calloc(M + 1, sizeof(double));
    double *pref_p = calloc((M + 1) * (M + 1), sizeof(double));
    double *b1 = calloc(M + 1, sizeof(double));
    double *b0_coeff = calloc(M + 1, sizeof(double));
    double B0_const = k * EminA;
    for (long d = 1; d < L && d <= M; d++)
        b0_coeff[d] = k * (2.0 / (double)(L * L)) * (double)(L - d);
    for (long y = 1; y <= M; y++) {
        long lt_low = A_low, lt_high = A_high;
        if (lt_high > y - 1) lt_high = y - 1;
        double sum_lt_a = arith_sum(lt_low, lt_high);
        if (lt_low <= lt_high) {
            for (long j = y - lt_high; j <= y - lt_low; j++)
                if (j >= 1 && j <= M) p0[y * (M + 1) + j] += invL;
        }
        long gt_low = A_low; if (gt_low < y + 1) gt_low = y + 1;
        long gt_high = A_high;
        long count_gt = gt_high - gt_low + 1; if (count_gt < 0) count_gt = 0;
        double base_gt = (double)count_gt * (double)y;
        if (gt_low <= gt_high) {
            for (long j = gt_low - y; j <= gt_high - y; j++)
                if (j >= 1 && j <= M) p0[y * (M + 1) + j] += invL;
        }
        double q = (sum_lt_a + base_gt) * invL;
        if (A_low <= y && y <= A_high) {
            q += (double)y * invL + B0_const * invL;
            for (long d = 1; d < L && d <= M; d++)
                p0[y * (M + 1) + d] += b0_coeff[d] * invL;
        }
        q0[y] = q;
    }
    prefix_fill(q0, pref_q, M);
    for (long j = 1; j <= M; j++) {
        double s = 0.0;
        for (long y = 1; y <= M; y++) {
            s += p0[y * (M + 1) + j];
            pref_p[j * (M + 1) + y] = s;
        }
    }
    double sum_w_q = 0.0;
    for (long d = 1; d < L && d <= M; d++) sum_w_q += (double)(L - d) * q0[d];
    double constB1 = k * EminA + k * (1.0 / (double)(L * L)) * sum_w_q;
    for (long j = 1; j <= M; j++) {
        double ss = 0.0;
        for (long d = 1; d < L && d <= M; d++) ss += (double)(L - d) * p0[d * (M + 1) + j];
        b1[j] = k * (1.0 / (double)(L * L)) * ss;
    }
    double *A_mat = calloc(M * (M + 1), sizeof(double));
    double *b_vec = calloc(M, sizeof(double));
    double *sol = calloc(M, sizeof(double));
    for (long t = 1; t <= M; t++) {
        long lt_low = A_low, lt_high = A_high;
        if (lt_high > t - 1) lt_high = t - 1;
        double sum_lt_a = arith_sum(lt_low, lt_high);
        long gt_low = A_low; if (gt_low < t + 1) gt_low = t + 1;
        long gt_high = A_high;
        long count_gt = gt_high - gt_low + 1; if (count_gt < 0) count_gt = 0;
        double base_gt = (double)count_gt * (double)t;
        double const1 = (sum_lt_a + base_gt) * invL;
        if (A_low <= t && t <= A_high) const1 += (double)t * invL;
        double Qgt = 0.0;
        if (gt_low <= gt_high) {
            long y_lo = gt_low - t, y_hi = gt_high - t;
            Qgt = pref_q[y_hi] - (y_lo > 1 ? pref_q[y_lo - 1] : 0.0);
            for (long j = 1; j <= M; j++) {
                double Pgt = pref_p[j * (M + 1) + y_hi] - (y_lo > 1 ? pref_p[j * (M + 1) + (y_lo - 1)] : 0.0);
                A_mat[(t - 1) * (M + 1) + (j - 1)] -= invL * Pgt;
            }
        }
        double rhs = const1 + invL * Qgt;
        if (A_low <= t && t <= A_high) rhs += invL * constB1;
        b_vec[t - 1] = rhs;
        A_mat[(t - 1) * (M + 1) + (t - 1)] += 1.0;
        if (A_low <= t && t <= A_high) {
            for (long j = 1; j <= M; j++)
                A_mat[(t - 1) * (M + 1) + (j - 1)] -= invL * b1[j];
        }
    }
    gauss_solve(A_mat, b_vec, M, sol);
    for (long j = 1; j <= M; j++) C1[j] = sol[j - 1];
    double ss2 = 0.0;
    for (long d = 1; d < L; d++) ss2 += (double)(L - d) * 2.0 * C1[d];
    *B0_out = k * (EminA + ss2 / (double)(L * L));
    free(q0); free(p0); free(pref_q); free(pref_p); free(b1); free(b0_coeff);
    free(A_mat); free(b_vec); free(sol);
}

static double expected_game(long m, long n) {
    long L = m - n + 1;
    long M = m + 5;
    double *C1 = calloc(M + 1, sizeof(double));
    double B0 = 0.0;
    int used = 0;
    if (L >= 25) used = compute_C1_iterative(n, m, C1, &B0);
    if (!used) compute_C1_direct(n, m, C1, &B0);
    double Emin0 = expected_min_consecutive(n, m);
    double s = B0 * (1.0 / (double)L);
    double denom = (double)(L * L);
    for (long d = 1; d < L; d++) s += (2.0 * (double)(L - d) / denom) * C1[d];
    free(C1);
    return Emin0 + s;
}

static double S(long k) {
    double total = 0.0;
    for (long m = 2; m <= k; m++)
        for (long n = 1; n < m; n++)
            total += expected_game(m, n);
    return total;
}

double pe589_answer(void) { return S(100); }
