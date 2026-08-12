// Project Euler 909: L-Expressions I
// F(n) = g(g(n^2*(n+1))) where g(x) = x*(x+1), mod 10^9
// seed = F(1) = 42, ans = F(seed)
// Print last 9 digits zero-padded.
#include <stdint.h>
#include <stdio.h>

#define MOD 1000000000LL

static long long F_mod(long long n, long long mod) {
    n %= mod;
    long long a = n * n % mod;
    a = a * ((n + 1) % mod) % mod;
    long long b = a * ((a + 1) % mod) % mod;
    long long c = b * ((b + 1) % mod) % mod;
    return c;
}

long long p909_native(void) {
    long long seed = F_mod(1, MOD);
    long long ans = F_mod(seed, MOD);
    return ans;
}
