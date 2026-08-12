// Project Euler 903: Total Permutation Powers
// Compute Q(n) mod 1e9+7 for n = 10^6.
// Uses Mobius sieve, harmonic numbers, and modular arithmetic.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MOD 1000000007LL

static long long mod_pow(long long base, long long exp, long long mod) {
    long long result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = result * base % mod;
        base = base * base % mod;
        exp >>= 1;
    }
    return result;
}

long long p903_native(void) {
    long long n = 1000000LL;

    // Modular inverses inv[1..n]
    long long *inv = (long long *)malloc((n + 1) * sizeof(long long));
    inv[1] = 1;
    for (long long i = 2; i <= n; i++) {
        inv[i] = MOD - (MOD / i) * inv[MOD % i] % MOD;
    }

    // Harmonic numbers H[k] = sum_{i=1..k} inv[i]
    long long *H = (long long *)malloc((n + 1) * sizeof(long long));
    long long s = 0;
    for (long long i = 1; i <= n; i++) {
        s += inv[i];
        if (s >= MOD) s -= MOD;
        H[i] = s;
    }

    // Mobius function via linear sieve
    signed char *mu = (signed char *)calloc(n + 1, 1);
    int *primes = (int *)malloc(n * sizeof(int));
    int prime_count = 0;
    char *is_comp = (char *)calloc(n + 1, 1);

    mu[1] = 1;
    for (int i = 2; i <= n; i++) {
        if (!is_comp[i]) {
            primes[prime_count++] = i;
            mu[i] = -1;
        }
        for (int j = 0; j < prime_count; j++) {
            int p = primes[j];
            int v = i * p;
            if ((long long)v > n) break;
            is_comp[v] = 1;
            if (i % p == 0) {
                mu[v] = 0;
                break;
            }
            mu[v] = -mu[i];
        }
    }

    // Compute F[s] = sum_{d|s} mu(d)/d * H[s/d - 1] by convolution over multiples
    long long *F = (long long *)calloc(n + 1, sizeof(long long));
    for (int d = 1; d <= n; d++) {
        if (mu[d] == 0) continue;
        long long c = (mu[d] == 1) ? inv[d] : (MOD - inv[d]);
        int m = 1;
        for (int s2 = d; s2 <= n; s2 += d) {
            F[s2] = (F[s2] + c * H[m - 1]) % MOD;
            m++;
        }
    }

    // S = sum_{s=2..n} H[floor(n/s)] * 2*F[s]/s
    long long S = 0;
    for (int s2 = 2; s2 <= n; s2++) {
        long long term = H[n / s2] * ((2 * F[s2]) % MOD) % MOD;
        term = term * inv[s2] % MOD;
        S = (S + term) % MOD;
    }

    // alpha = (n - H[n] + S) / (n * (n-1))
    long long num_alpha = ((n % MOD - H[n] + S) % MOD + MOD) % MOD;
    long long alpha = num_alpha * inv[n] % MOD * inv[n - 1] % MOD;

    // beta = H[n/2] / (2*n*(n-1))
    long long denom_beta = (2 * n % MOD) * (n - 1) % MOD;
    long long beta = H[n / 2] * mod_pow(denom_beta, MOD - 2, MOD) % MOD;

    // p = H[n] / n
    long long p = H[n] * inv[n] % MOD;
    // q = (1 - p) / (n-1)
    long long q = ((1 - p) % MOD + MOD) % MOD * inv[n - 1] % MOD;

    long long inv_n2 = inv[n - 2];
    long long inv_n3 = inv[n - 3];

    // a = (p - alpha) / (n-2)
    long long a = ((p - alpha) % MOD + MOD) % MOD * inv_n2 % MOD;
    // b = (q - beta) / (n-2)
    long long b = ((q - beta) % MOD + MOD) % MOD * inv_n2 % MOD;
    // eta = (q - a - b) / (n-3)
    long long eta = ((q - a - b) % MOD + MOD) % MOD * inv_n3 % MOD;

    // C0 = beta + (n-3)*b + (n-1)*a + eta * (n-2)(n-3)/2
    long long bconst = ((n - 2) * (n - 3) / 2) % MOD;
    long long C0 = (beta + (n - 3) % MOD * b % MOD + (n - 1) % MOD * a % MOD + eta * bconst % MOD) % MOD;
    C0 = (C0 % MOD + MOD) % MOD;
    long long slope = ((b - a) % MOD + MOD) % MOD;

    long long inv2 = (MOD + 1) / 2;

    // S1 = sum_{m=1..n-1} m! * m
    // S2 = sum_{m=1..n-1} m! * m*(m+1)/2
    long long fact = 1, S1 = 0, S2 = 0;
    for (int m = 1; m < n; m++) {
        fact = fact * m % MOD;  // fact = m!
        S1 = (S1 + fact * m) % MOD;
        S2 = (S2 + fact * m % MOD * (m + 1) % MOD * inv2) % MOD;
    }

    long long fact_n = fact * n % MOD;  // n!
    long long E_rank = (1 + C0 * S1 % MOD + slope * S2 % MOD) % MOD;
    E_rank = (E_rank % MOD + MOD) % MOD;

    long long ans = fact_n * fact_n % MOD * E_rank % MOD;

    free(inv); free(H); free(mu); free(primes); free(is_comp); free(F);
    return ans;
}
