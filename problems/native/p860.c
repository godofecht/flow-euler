// Project Euler 860: Gold and Silver Coin Game
// Count fair arrangements of n stacks mod 989898989.
#include <stdint.h>
#include <stdlib.h>
typedef long long i64;

enum { MOD = 989898989 };
#define N 9898

static i64 modinv(i64 a, i64 mod) {
    a %= mod;
    if (a < 0) a += mod;
    i64 t0 = 0, t1 = 1, r0 = mod, r1 = a;
    while (r1) {
        i64 q = r0 / r1;
        i64 tmp;
        tmp = r0 - q * r1; r0 = r1; r1 = tmp;
        tmp = t0 - q * t1; t0 = t1; t1 = tmp;
    }
    return t0 % mod;
}

long long p860_native(void) {
    i64 *fac = malloc((N + 1) * sizeof(i64));
    i64 *inv_fac = malloc((N + 1) * sizeof(i64));
    fac[0] = 1;
    for (int i = 1; i <= N; i++)
        fac[i] = fac[i - 1] * i % MOD;
    inv_fac[N] = modinv(fac[N], MOD);
    for (int i = N; i >= 1; i--)
        inv_fac[i - 1] = inv_fac[i] * i % MOD;

    i64 fac_n = fac[N];
    i64 total = 0;

    for (int t = 0; t <= N; t++) {
        int s = N - t;
        if (s & 1) continue;
        i64 choose_nt = fac_n * inv_fac[t] % MOD * inv_fac[s] % MOD;
        int high = (t < s / 4) ? t : (s / 4);
        int parity = t & 1;
        i64 inner = 0;
        i64 fac_t = fac[t];
        i64 fac_s = fac[s];

        for (int m = parity; m <= high; m += 2) {
            int k = (t + m) / 2;
            int r = (s - 4 * m) / 2;
            i64 ct = fac_t * inv_fac[k] % MOD * inv_fac[t - k] % MOD;
            i64 cs = fac_s * inv_fac[r] % MOD * inv_fac[s - r] % MOD;
            i64 term = ct * cs % MOD;
            if (m == 0)
                inner += term;
            else
                inner += term * 2;
        }
        total = (total + choose_nt * (inner % MOD)) % MOD;
    }

    free(fac);
    free(inv_fac);
    return total;
}
