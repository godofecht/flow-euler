/* Project Euler 959: Asymmetric Random Walk.
 *
 * A frog starts at 0. Each step jumps a units left or b units right with
 * probability 1/2 each. f(a,b) = lim c_n/n where c_n is the expected number
 * of distinct positions visited in n steps.
 *
 * For a != b the walk is transient and f(a,b) = 1/G where G is the expected
 * number of visits to the origin:
 *   G = sum_{k>=0} C((a+b)k, ak) / 2^{(a+b)k}
 *
 * Each term is computed in log space via lgamma to avoid overflow/underflow
 * in the intermediate factorials.
 */

#include <math.h>

static long gcd_long(long x, long y) {
    while (y != 0) {
        long t = x % y;
        x = y;
        y = t;
    }
    return x < 0 ? -x : x;
}

static double f_ab(long a, long b) {
    long g = gcd_long(a, b);
    a /= g;
    b /= g;

    if (a == b) return 0.0;
    if (a > b) { long t = a; a = b; b = t; }

    long m = a + b;
    double log2 = log(2.0);

    /* t_k = C(mk, ak) / 2^{mk}
     *     = exp(lgamma(mk+1) - lgamma(ak+1) - lgamma(bk+1) - mk*log 2)
     */
    double S = 0.0;
    long k = 0;
    while (1) {
        long mk = m * k;
        long ak = a * k;
        long bk = b * k;

        double log_t = lgamma((double)(mk + 1))
                     - lgamma((double)(ak + 1))
                     - lgamma((double)(bk + 1))
                     - (double)mk * log2;

        if (log_t < -700.0) break;  /* term underflows double; tail is negligible */

        double t = exp(log_t);
        S += t;

        if (t < 1e-17 * S) break;
        k++;
        if (k > 200000) break;
    }

    return 1.0 / S;
}

double p959_native(void) {
    return f_ab(89, 97);
}
