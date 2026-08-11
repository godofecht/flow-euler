/*
 * Project Euler 756: Approximating a Sum
 *
 * Computes E(Delta | phi(k), n, m) exactly (no simulation) for
 * n = 12345678, m = 12345, where f(k) = phi(k) (Euler's totient).
 *
 * S  = sum_{k=1..n} f(k)
 * S* = sum_{i=1..m} f(X_i) (X_i - X_{i-1})   for a random increasing m-tuple
 * Delta = S - S*
 *
 * E(Delta) = sum_{k=1..n-m} f(k) * C(n-k, m) / C(n, m)
 * with weights updated by recurrence w_{k+1} = w_k * (n-k-m)/(n-k).
 *
 * A rigorous tail bound truncates the sum so only the necessary initial
 * part is computed, keeping the totient sieve small.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Euler's linear sieve: phi[k] = totient(k) for 0 <= k <= n. */
static int *totients_up_to(long n) {
    int *phi = (int *)calloc((size_t)n + 1, sizeof(int));
    if (!phi) {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }
    if (n >= 1) phi[1] = 1;

    /* primes list */
    long cap = 0;
    long count = 0;
    int *primes = NULL;

    for (long i = 2; i <= n; i++) {
        if (phi[i] == 0) { /* i is prime */
            if (count == cap) {
                cap = cap ? cap * 2 : 256;
                primes = (int *)realloc(primes, (size_t)cap * sizeof(int));
                if (!primes) { fprintf(stderr, "oom\n"); exit(1); }
            }
            primes[count++] = (int)i;
            phi[i] = (int)(i - 1);
        }
        for (long j = 0; j < count; j++) {
            long p = primes[j];
            long ip = i * p;
            if (ip > n) break;
            if (i % p == 0) {
                phi[ip] = phi[i] * (int)p;
                break;
            } else {
                phi[ip] = phi[i] * (int)(p - 1);
            }
        }
    }
    free(primes);
    return phi;
}

/* Find a safe truncation index K (<= n-m) such that the remaining tail
 * contribution is provably < eps. Uses phi(k) <= k <= n and non-increasing
 * weights: sum_{k>K} phi(k) w_k <= n * (n-m-K) * w_{K+1}. */
static long cutoff_index(long n, long m, double eps) {
    long limit = n - m;
    if (limit <= 0) return 0;

    double w = (double)(n - m) / (double)n; /* w_1 */
    for (long k = 1; k <= limit; k++) {
        long remaining = limit - k;
        long nk = n - k;
        double w_next;
        if (nk <= m) {
            w_next = 0.0;
        } else {
            w_next = w * (double)(nk - m) / (double)nk;
        }
        if ((double)n * (double)remaining * w_next < eps) {
            return k;
        }
        w = w_next;
    }
    return limit;
}

static double expected_error_for_phi(long n, long m) {
    long limit = n - m;
    if (limit <= 0) return 0.0;

    long upto = cutoff_index(n, m, 5e-8);
    if (upto > limit) upto = limit;

    int *phi = totients_up_to(upto);

    double w = (double)(n - m) / (double)n; /* w_1 */
    double total = 0.0;
    double c = 0.0; /* Kahan compensation */

    for (long k = 1; k <= upto; k++) {
        double x = (double)phi[k] * w;
        /* Kahan compensated summation */
        double y = x - c;
        double t = total + y;
        c = (t - total) - y;
        total = t;

        long nk = n - k;
        if (nk <= m) break;
        w *= (double)(nk - m) / (double)nk;
    }

    free(phi);
    return total;
}

double p756_native(void) {
    long n = 12345678L;
    long m = 12345L;
    return expected_error_for_phi(n, m);
}
