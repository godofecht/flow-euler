#include <stdint.h>
#include <stdlib.h>
typedef long long i64;
enum { MOD = 1000000007LL };
static i64 modpow(i64 b, i64 e) {
    i64 r = 1; b %= MOD;
    while (e > 0) { if (e & 1) r = r * b % MOD; b = b * b % MOD; e >>= 1; }
    return r;
}
long long pe_solve(void) {
    int N = 2022;
    i64 *fac = calloc(N + 1, sizeof(i64));
    fac[0] = 1;
    for (int x = 1; x <= N; x++) fac[x] = fac[x - 1] * x % MOD;
    i64 *invfac = calloc(N + 1, sizeof(i64));
    invfac[N] = modpow(fac[N], MOD - 2);
    for (int x = N; x >= 1; x--) invfac[x - 1] = invfac[x] * x % MOD;
    i64 total = 0;
    for (int n = 1; n <= N; n++) {
        for (int m = n / 2 + 1; m <= n; m++) {
            i64 c = modpow(9, n - m + 1) * fac[n] % MOD * invfac[n - m] % MOD * invfac[m] % MOD;
            total = (total + c) % MOD;
        }
    }
    free(fac); free(invfac);
    return total;
}
