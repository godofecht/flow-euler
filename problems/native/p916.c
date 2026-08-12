// Project Euler 916: Restricted Permutations
// P(n) = Catalan(n)^2 * (1 + (3n/(n+2))^2) mod 1e9+7
// n = 10^8
#include <stdint.h>
#include <stdio.h>

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

long long p916_native(void) {
    long long n = 100000000LL;
    // Compute n! mod MOD and (2n)! mod MOD
    long long fac_n = 1, acc = 1;
    for (long long i = 1; i <= 2 * n; i++) {
        acc = acc * (i % MOD) % MOD;
        if (i == n) fac_n = acc;
    }
    long long fac_2n = acc;

    long long inv_fac_n = mod_pow(fac_n, MOD - 2, MOD);
    long long binom_2n_n = fac_2n * inv_fac_n % MOD * inv_fac_n % MOD;
    long long catalan = binom_2n_n * mod_pow(n + 1, MOD - 2, MOD) % MOD;

    long long t = 3 * n % MOD * mod_pow(n + 2, MOD - 2, MOD) % MOD;

    long long ans = catalan * catalan % MOD * ((1 + t * t % MOD) % MOD) % MOD;
    return ans;
}
