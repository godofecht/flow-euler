/*
 * Project Euler 779
 * Prime factor and exponent.
 *
 * sum_{K>=1} bar(f_K) =
 *   sum_{p prime} prod_{q<p}(1 - 1/q) * 1 / (p * (p - 1)^2)
 *
 * Reference: /tmp/pes_ref/solvers/779.py (sum_f).
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SIEVE_LIMIT 1000000

double p779_native(void) {
    /* Sieve of Eratosthenes up to SIEVE_LIMIT. */
    char *is_prime = (char *)malloc((size_t)(SIEVE_LIMIT + 1));
    if (!is_prime) {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }
    memset(is_prime, 1, (size_t)(SIEVE_LIMIT + 1));
    is_prime[0] = 0;
    is_prime[1] = 0;
    long limit = (long)floor(sqrt((double)SIEVE_LIMIT));
    for (long i = 2; i <= limit; i++) {
        if (is_prime[i]) {
            for (long j = i * i; j <= SIEVE_LIMIT; j += i) {
                is_prime[j] = 0;
            }
        }
    }

    double total = 0.0;
    double curr = 1.0; /* running product of (p-1)/p over primes seen */
    for (long p = 2; p <= SIEVE_LIMIT; p++) {
        if (!is_prime[p]) continue;
        double dp = (double)p;
        double term = curr * (1.0 / (dp * (dp - 1.0) * (dp - 1.0)));
        total += term;
        curr *= (dp - 1.0) / dp;
    }

    free(is_prime);

    /* Round to 12 decimal places like the reference. */
    return round(total * 1e12) / 1e12;
}
