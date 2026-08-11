#include <stdlib.h>

/* Project Euler 819 - Iterative Sampling.
 *
 * Starting from the n-tuple (1,2,...,n), each step creates a new n-tuple by
 * independently sampling each position from the multiset of values in the
 * previous tuple. E(n) is the expected number of steps until all entries are
 * equal.
 *
 * The process is equivalent to repeated composition of random functions
 * f:[n]->[n]. The only state needed is m = size of the image of the composed
 * function (number of distinct ancestral labels still present). Given m, the
 * next image size is the number of occupied bins when throwing m balls
 * uniformly into n bins.
 *
 * The chain is monotone non-increasing in m, so the expected hitting time to
 * m=1 is computed by a triangular dynamic program.
 */

double p819_native(void) {
    int n = 1000;
    if (n <= 1) return 0.0;

    /* T[m] = expected remaining steps to reach 1 distinct label, starting from m */
    double *T = (double *)calloc((size_t)(n + 1), sizeof(double));
    T[1] = 0.0;

    /* p[k] = P(number of occupied bins = k) for the current m (updated in place) */
    double *p = (double *)calloc((size_t)(n + 1), sizeof(double));
    p[0] = 1.0;
    double inv_n = 1.0 / (double)n;

    for (int m = 1; m <= n; m++) {
        /* Update occupancy distribution from m-1 balls to m balls:
         * new_p[k] = old_p[k] * (k/n) + old_p[k-1] * ((n-k+1)/n)
         * Do it in-place descending so old_p[k-1] is still available. */
        for (int k = m; k >= 1; k--) {
            p[k] = p[k] * ((double)k * inv_n) + p[k - 1] * ((double)(n - k + 1) * inv_n);
        }
        p[0] = 0.0;

        if (m == 1) continue;

        double stay = p[m]; /* probability of no collisions -> state stays at m */
        double acc = 0.0;
        for (int k = 1; k < m; k++) {
            acc += p[k] * T[k];
        }

        /* T[m] = 1 + stay*T[m] + sum_{k<m} p[k]*T[k] */
        T[m] = (1.0 + acc) / (1.0 - stay);
    }

    double result = T[n];
    free(T);
    free(p);
    return result;
}
