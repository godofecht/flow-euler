#include <math.h>
#include <stdio.h>

/* Project Euler 825 - Chasing Game.
 *
 * Computes T(10^14) rounded to 8 decimal places.
 *
 * The algorithm mirrors the Python reference:
 *   c = (3 - sqrt(3)) / 6
 *   K = correction constant (sum of S(n) - 1/(n-1+c) for n=2..60)
 *   T(N) = shifted_harmonic_sum(N, c) + K
 *      = (digamma(N+c) - digamma(1+c)) + K
 *
 * S_fast_float(n) uses a closed form with q = -2+sqrt(3) and a 4x4 linear system.
 */

static double digamma(double x) {
    double res = 0.0;
    while (x < 8.0) {
        res -= 1.0 / x;
        x += 1.0;
    }
    double inv = 1.0 / x;
    double inv2 = inv * inv;
    res += log(x) - 0.5 * inv;
    double t = inv2;
    res -= t / 12.0;
    t *= inv2;
    res += t / 120.0;
    t *= inv2;
    res -= t / 252.0;
    t *= inv2;
    res += t / 240.0;
    t *= inv2;
    res -= t / 132.0;
    t *= inv2;
    res += t * (691.0 / 32760.0);
    return res;
}

static void gauss_solve_4x4(double mat[4][5], double sol[4]) {
    for (int c = 0; c < 4; c++) {
        int piv = c;
        for (int r = c + 1; r < 4; r++) {
            if (fabs(mat[r][c]) > fabs(mat[piv][c])) piv = r;
        }
        if (piv != c) {
            for (int j = 0; j < 5; j++) {
                double tmp = mat[c][j];
                mat[c][j] = mat[piv][j];
                mat[piv][j] = tmp;
            }
        }
        double pv = mat[c][c];
        for (int j = c; j < 5; j++) mat[c][j] /= pv;
        for (int r = 0; r < 4; r++) {
            if (r == c) continue;
            double f = mat[r][c];
            if (f == 0.0) continue;
            for (int j = c; j < 5; j++) mat[r][j] -= f * mat[c][j];
        }
    }
    for (int i = 0; i < 4; i++) sol[i] = mat[i][4];
}

static double S_fast_float(int n) {
    if (n < 2) return 0.0;
    if (n == 2) {
        /* S_fraction(2) = 7/11, so 2*(7/11)-1 = 3/11 */
        return 7.0 / 11.0;
    }
    int L = 2 * n;
    double q = -2.0 + sqrt(3.0);

    /* row(y) = [1, y, q^y, q^(L-y)] */
    double row[6][4]; /* for y = 2,3,4,5,L-2,L-1 */
    int ys[6] = {2, 3, 4, 5, L - 2, L - 1};
    for (int i = 0; i < 6; i++) {
        int y = ys[i];
        row[i][0] = 1.0;
        row[i][1] = (double)y;
        row[i][2] = pow(q, (double)y);
        row[i][3] = pow(q, (double)(L - y));
    }

    /* indices: row(2)=0, row(3)=1, row(4)=2, row(5)=3, row(L-2)=4, row(L-1)=5 */
    double mat[4][5];
    /* Eq1: g(2) + (1/3) g(L-1) = 1 */
    for (int j = 0; j < 4; j++) mat[0][j] = row[0][j] + (1.0/3.0) * row[5][j];
    mat[0][4] = 1.0;
    /* Eq2: g(3) + (1/3) g(L-2) + (1/3) g(L-1) = 1 */
    for (int j = 0; j < 4; j++) mat[1][j] = row[1][j] + (1.0/3.0) * row[4][j] + (1.0/3.0) * row[5][j];
    mat[1][4] = 1.0;
    /* Eq3: g(L-1) + (1/3)g(2) + (1/3)g(3) + (1/3)g(4) = 1 */
    for (int j = 0; j < 4; j++) mat[2][j] = row[5][j] + (1.0/3.0) * row[0][j] + (1.0/3.0) * row[1][j] + (1.0/3.0) * row[2][j];
    mat[2][4] = 1.0;
    /* Eq4: g(L-2) + (1/3)g(3) + (1/3)g(4) + (1/3)g(5) = 1 */
    for (int j = 0; j < 4; j++) mat[3][j] = row[4][j] + (1.0/3.0) * row[1][j] + (1.0/3.0) * row[2][j] + (1.0/3.0) * row[3][j];
    mat[3][4] = 1.0;

    double sol[4];
    gauss_solve_4x4(mat, sol);
    double A = sol[0], B = sol[1], C = sol[2], E = sol[3];

    double gn = A + B * (double)n + C * pow(q, (double)n) + E * pow(q, (double)(L - n));
    return 2.0 * gn - 1.0;
}

static double correction_constant(void) {
    double c = (3.0 - sqrt(3.0)) / 6.0;
    double acc = 0.0;
    for (int n = 2; n <= 60; n++) {
        acc += S_fast_float(n) - 1.0 / ((double)(n - 1) + c);
    }
    return acc;
}

double p825_native(void) {
    double c = (3.0 - sqrt(3.0)) / 6.0;
    double K = correction_constant();
    double N = 1e14;
    double shs = digamma(N + c) - digamma(1.0 + c);
    double ans = shs + K;
    return ans;
}
