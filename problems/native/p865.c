// Project Euler 865: Triplicate Numbers
// T(10^4) mod 998244353.
#include <stdint.h>
#include <stdlib.h>
typedef long long i64;

enum { MOD = 998244353 };
#define K 10
#define MMAX (10000 / 3)  // 3333

long long p865_native(void) {
    int mmax = MMAX;
    i64 *dp0 = calloc(mmax + 1, sizeof(i64));
    i64 *dp1 = calloc(mmax + 1, sizeof(i64));
    i64 *f = calloc(mmax + 1, sizeof(i64));
    i64 *mul = calloc(mmax + 1, sizeof(i64));
    i64 *pref = calloc(mmax + 1, sizeof(i64));

    i64 k = K, km1 = K - 1, twok = 2 * K, twokm1 = 2 * (K - 1);

    dp0[1] = k;
    dp1[1] = km1;
    f[1] = km1;
    pref[1] = dp1[1];

    for (int m = 2; m <= mmax; m++) {
        i64 f_m = km1 * f[m - 1] % MOD;
        i64 dp0_m = k * dp0[m - 1] % MOD;
        i64 dp1_m = k * dp1[m - 1] % MOD;

        for (int s = 2; s <= m; s++) {
            int p = m - s;
            i64 x = f[s - 1];

            if (p == 0) {
                f_m = (f_m + twokm1 * x) % MOD;
                dp0_m = (dp0_m + twok * x) % MOD;
                dp1_m = (dp1_m + twokm1 * x) % MOD;
                if (s >= 3) {
                    i64 y = mul[s - 1];
                    f_m = (f_m + km1 * y) % MOD;
                    dp0_m = (dp0_m + k * y) % MOD;
                    dp1_m = (dp1_m + km1 * y) % MOD;
                }
            } else {
                i64 fp = f[p];
                i64 d0p = dp0[p];
                i64 d1p = dp1[p];
                f_m = (f_m + twokm1 * x % MOD * fp) % MOD;
                dp0_m = (dp0_m + twok * x % MOD * d0p) % MOD;
                dp1_m = (dp1_m + twok * x % MOD * d1p) % MOD;
                if (s >= 3) {
                    i64 y = mul[s - 1];
                    f_m = (f_m + km1 * y % MOD * fp) % MOD;
                    dp0_m = (dp0_m + k * y % MOD * d0p) % MOD;
                    dp1_m = (dp1_m + k * y % MOD * d1p) % MOD;
                }
            }

            if ((s & 63) == 0) {
                f_m %= MOD;
                dp0_m %= MOD;
                dp1_m %= MOD;
            }
        }

        f_m %= MOD;
        dp0_m %= MOD;
        dp1_m %= MOD;

        f[m] = f_m;
        dp0[m] = dp0_m;
        dp1[m] = dp1_m;

        i64 acc = 0;
        for (int a = 1; a < m; a++) {
            acc = (acc + f[a] * f[m - a] % MOD) % MOD;
        }
        mul[m] = acc;

        pref[m] = (pref[m - 1] + dp1[m]) % MOD;
    }

    i64 ans = pref[mmax] % MOD;
    free(dp0); free(dp1); free(f); free(mul); free(pref);
    return ans;
}
