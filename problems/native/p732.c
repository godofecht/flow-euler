/* Project Euler 732 - Standing on the Shoulders of Trolls
   Compute Q(1000): maximum total IQ of trolls that can escape.
   Generates trolls from r[i] = (5^i mod 1e9+7) mod 101 + 50, then solves
   a 0/1 knapsack-over-time scheduling problem (earliest-deadline-first DP).
*/
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MOD 1000000007LL

/* Integer square root for non-negative long long. */
static long long isqrt_ll(long long x) {
    if (x < 0) return 0;
    if (x < 2) return x;
    long long x0 = x;
    long long x1 = (x0 + 1) / 2;
    while (x1 < x0) {
        x0 = x1;
        x1 = (x0 + x / x0) / 2;
    }
    return x0;
}

/* ceil(s / sqrt(2)) exactly: smallest y with 2*y^2 >= s^2. */
static long long ceil_div_sqrt2(long long s) {
    long long a = s * s;
    return isqrt_ll((a - 1) / 2) + 1;
}

typedef struct {
    long long d;   /* deadline */
    long long p;   /* processing time (h) */
    long long profit; /* q */
} Job;

static int cmp_job(const void *a, const void *b) {
    const Job *ja = (const Job *)a;
    const Job *jb = (const Job *)b;
    if (ja->d < jb->d) return -1;
    if (ja->d > jb->d) return 1;
    return 0;
}

long long p732_native(void) {
    int n = 1000;
    int total_r = 3 * n;

    /* Generate r[0..3n-1]. */
    long long *r = (long long *)malloc(sizeof(long long) * (size_t)total_r);
    long long p = 1; /* 5^0 mod MOD */
    for (int i = 0; i < total_r; i++) {
        r[i] = (p % 101) + 50;
        p = (p * 5) % MOD;
    }

    long long total_h = 0;
    for (int k = 0; k < n; k++) {
        total_h += r[3 * k];
    }

    long long y = ceil_div_sqrt2(total_h);
    long long base = total_h - y;

    /* Build jobs. */
    Job *jobs = (Job *)malloc(sizeof(Job) * (size_t)n);
    int njobs = 0;
    long long max_d = 0;
    for (int k = 0; k < n; k++) {
        long long h = r[3 * k];
        long long l = r[3 * k + 1];
        long long q = r[3 * k + 2];
        long long d = base + l + h;
        if (h <= d) {
            jobs[njobs].d = d;
            jobs[njobs].p = h;
            jobs[njobs].profit = q;
            njobs++;
            if (d > max_d) max_d = d;
        }
    }

    qsort(jobs, (size_t)njobs, sizeof(Job), cmp_job);

    /* 0/1 knapsack DP over time. dp[t] = max profit with total time exactly t. */
    long long sz = max_d + 1;
    long long *dp = (long long *)malloc(sizeof(long long) * (size_t)sz);
    for (long long i = 0; i < sz; i++) dp[i] = -1;
    dp[0] = 0;

    for (int j = 0; j < njobs; j++) {
        long long d = jobs[j].d;
        long long pp = jobs[j].p;
        long long profit = jobs[j].profit;
        for (long long t = d; t >= pp; t--) {
            if (dp[t - pp] != -1) {
                long long cand = dp[t - pp] + profit;
                if (cand > dp[t]) dp[t] = cand;
            }
        }
    }

    long long best = 0;
    for (long long i = 0; i < sz; i++) {
        if (dp[i] > best) best = dp[i];
    }

    free(dp);
    free(jobs);
    free(r);
    return best;
}
