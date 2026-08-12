// Project Euler 874: Maximal Prime Score
// M(k, n): choose a_i in [0, k) with sum(a_i) divisible by k, maximize sum(p(a_i)).
// Strategy: start from all a_i = k-1, find minimal reduction r = (-n) % k,
// then minimize prime-score loss via unbounded knapsack.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int *first_n_primes(int n) {
    if (n <= 0) return NULL;
    int limit;
    if (n < 6) limit = 15;
    else {
        double nn = (double)n;
        limit = (int)(nn * (log(nn) + log(log(nn))) + 10);
    }

    while (1) {
        char *sieve = malloc(limit + 1);
        memset(sieve, 1, limit + 1);
        sieve[0] = sieve[1] = 0;
        for (int i = 2; (long long)i * i <= limit; i++) {
            if (sieve[i]) {
                for (int j = i * i; j <= limit; j += i)
                    sieve[j] = 0;
            }
        }
        int count = 0;
        for (int i = 2; i <= limit; i++)
            if (sieve[i]) count++;
        if (count >= n) {
            int *primes = malloc(n * sizeof(int));
            int idx = 0;
            for (int i = 2; i <= limit && idx < n; i++)
                if (sieve[i]) primes[idx++] = i;
            free(sieve);
            return primes;
        }
        free(sieve);
        limit *= 2;
    }
}

// Unbounded knapsack: minimal total loss to achieve exact reduction sum.
static long long min_prime_loss(const int *primes, int k, int reduction) {
    if (reduction == 0) return 0;
    int m = k - 1; // index of top prime
    int top = primes[m];

    // loss[d] = top - primes[m - d] for d = 1..reduction
    long long *loss = malloc((reduction + 1) * sizeof(long long));
    for (int d = 1; d <= reduction; d++)
        loss[d] = top - primes[m - d];

    // dp[s] = minimal loss to achieve sum s
    long long INF = (long long)1e18;
    long long *dp = malloc((reduction + 1) * sizeof(long long));
    dp[0] = 0;
    for (int s = 1; s <= reduction; s++)
        dp[s] = INF;

    for (int d = 1; d <= reduction; d++) {
        long long c = loss[d];
        for (int s = d; s <= reduction; s++) {
            long long cand = dp[s - d] + c;
            if (cand < dp[s])
                dp[s] = cand;
        }
    }

    long long result = dp[reduction];
    free(loss);
    free(dp);
    return result;
}

static long long maximal_prime_score(int k, long long n, const int *primes) {
    if (k == 1) return n * 2;
    int top = primes[k - 1];
    long long r = ((-n) % k + k) % k; // (-n) mod k, non-negative
    long long loss = min_prime_loss(primes, k, (int)r);
    return n * top - loss;
}

long long p874_native(void) {
    int k = 7000;
    int *primes_7001 = first_n_primes(k + 1);
    long long n = primes_7001[k]; // p(7000)
    long long ans = maximal_prime_score(k, n, primes_7001);
    free(primes_7001);
    return ans;
}
