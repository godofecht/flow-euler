// Project Euler 770 - Delphi Flip
//
// g(X) is the smallest n such that player A can guarantee at least X grams.
// Closed form: F(n) = 2 / (1 + p_n), where p_n = C(2n,n) / 4^n.
// The condition F(n) >= X is equivalent to p_n <= (2 - X) / X.
//
// For modest n we compare in log-space via lgamma (exact to double precision).
// For huge n we use the Stirling-series expansion of ln p_n, which avoids the
// catastrophic cancellation between two lgamma values of ~1e9 magnitude:
//
//     ln p_n = -1/2 ln(pi n) - 1/(8n) + 1/(192 n^3) - 1/(640 n^5) + O(n^-7)

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double ln_p_stirling(double n) {
    double inv = 1.0 / n;
    double inv2 = inv * inv;
    double inv3 = inv2 * inv;
    double inv5 = inv3 * inv2;
    return -0.5 * log(M_PI * n)
         - 0.125 * inv
         + (1.0 / 192.0) * inv3
         - (1.0 / 640.0) * inv5;
}

// ln(C(2n,n)/4^n) via lgamma.  Fine for small n where no cancellation occurs.
static double ln_p_lgamma(long n) {
    double nf = (double)n;
    return lgamma(2.0 * nf + 1.0)
         - 2.0 * lgamma(nf + 1.0)
         - 2.0 * nf * log(2.0);
}

static long g_for_fraction(long x_num, long x_den) {
    // r = (2 - X) / X = (2*x_den - x_num) / x_num
    long r_num = 2 * x_den - x_num;
    long r_den = x_num;

    double r = (double)r_num / (double)r_den;
    double ln_r = log((double)r_num) - log((double)r_den);

    // p_n ~ 1/sqrt(pi n)  =>  n ~ 1/(pi r^2)
    long n_est = (long)(1.0 / (M_PI * r * r));

    if (n_est < 20000) {
        long n = 0;
        while (ln_p_lgamma(n) > ln_r) n++;
        return n;
    }

    // Huge n: walk in log-space using the Stirling expansion.
    long n = n_est - 10;
    if (n < 1) n = 1;
    while (ln_p_stirling((double)n) > ln_r) n++;
    while (n > 1 && ln_p_stirling((double)(n - 1)) <= ln_r) n--;
    return n;
}

long long p770_native(void) {
    return (long long)g_for_fraction(19999, 10000);
}
