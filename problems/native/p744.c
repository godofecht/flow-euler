/*
 * Project Euler 744: What? Where? When?
 *
 * The order in which envelopes are opened is a uniformly random permutation.
 * Let Q be the number of questions needed for either side to reach n points
 * in an i.i.d. Bernoulli(p) sequence (ignoring the RED card).
 *
 * The RED card position K is uniform on {1,2,...,2n+1} and independent of
 * answers.  The game ends normally iff K > Q, hence:
 *
 *     f(n,p) = P(K > Q) = E[(2n+1 - Q)/(2n+1)] = 1 - E[Q]/(2n+1).
 *
 * For very large n with noticeable bias, the underdog winning is
 * astronomically unlikely, so E[Q] is essentially the negative binomial
 * mean n / max(p, 1-p).
 *
 * Reference: /tmp/pes_ref/solvers/744.py
 */

#include <math.h>

double p744_native(void) {
    /* n = 10^11, p = 0.4999 */
    long long n = 100000000000LL;
    double p = 0.4999;
    double q = 1.0 - p;
    double a = (p >= q) ? p : q;

    /* E[Q] ~= n / max(p, 1-p) for huge n. */
    double eq = (double)n / a;
    double denom = (double)(2 * n + 1);
    return 1.0 - eq / denom;
}
