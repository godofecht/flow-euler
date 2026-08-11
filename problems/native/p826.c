
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Project Euler 826 - Birds on a Wire.
 *
 * Conjecture: for a prime p, the expected painted fraction is
 *   F(p) = (7p + 15) / (18(p + 1)).
 *
 * The answer is the average of F(p) over all primes p with 3 <= p <= 10^6,
 * rounded to 10 decimal places.
 */

static char *sieve_to(int limit, int *count_out) {
    char *is_prime = (char *)malloc((size_t)(limit + 1));
    memset(is_prime, 1, (size_t)(limit + 1));
    is_prime[0] = 0;
    is_prime[1] = 0;
    for (int i = 2; (long)i * i <= limit; i++) {
        if (is_prime[i]) {
            for (long j = (long)i * i; j <= limit; j += i) {
                is_prime[j] = 0;
            }
        }
    }
    int count = 0;
    for (int i = 2; i <= limit; i++) if (is_prime[i]) count++;
    *count_out = count;
    return is_prime;
}

static double F(long p) {
    double num = (double)(7 * p + 15);
    double den = (double)(18 * (p + 1));
    return num / den;
}

double p826_native(void) {
    int limit = 1000000;
    int prime_count = 0;
    char *is_prime = sieve_to(limit, &prime_count);

    /* primes[0] is 2; skip it, average over the rest. */
    double total = 0.0;
    for (int p = 3; p <= limit; p++) {
        if (is_prime[p]) total += F(p);
    }
    free(is_prime);

    double denom = (double)(prime_count - 1);
    double avg = total / denom;

    /* Round to 10 decimal places. */
    double scale = 10000000000.0;
    double rounded = round(avg * scale) / scale;
    return rounded;
}
