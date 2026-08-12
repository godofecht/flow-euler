// Project Euler 977: Count functions f:[n]->[n] with f^{(x)}(y)=f^{(y)}(x), mod 1e9+7.
// Port of the Python reference solver to C.
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;
typedef __int128 i128;

#define MOD 1000000007LL

static i64 mod_pow(i64 a, i64 e) {
    i64 r = 1;
    a %= MOD;
    if (a < 0) a += MOD;
    while (e > 0) {
        if (e & 1) r = (i64)((i128)r * a % MOD);
        a = (i64)((i128)a * a % MOD);
        e >>= 1;
    }
    return r;
}

long long p977_native(void) {
    i64 n = 1000000;
    if (n == 1) return 1;

    i64 total = 0;

    /* L = 1 closed form */
    i64 m = n - 2;
    i128 sum_q = (i128)m * (m + 1) * (2 * m + 1) / 6 + (i128)m * (m + 1) / 2;
    total = (i64)((sum_q % MOD + n) % MOD);

    /* Preallocate powA buffer: max size is for L=2, max_a ~ n/2 + 2 */
    i64 max_buf = n / 2 + 3;
    i64 *powA = (i64 *)malloc((size_t)max_buf * sizeof(i64));

    for (i64 L = 2; L <= n; L++) {
        i64 R = n - L;

        if (R >= 1) {
            i64 q_full = (R - 1) / L;
            i64 max_a = q_full + 2;

            /* Precompute a^L mod MOD for a in [1..max_a] */
            if (L == 2) {
                for (i64 a = 1; a <= max_a; a++)
                    powA[a] = (a * a) % MOD;
            } else if (L == 3) {
                for (i64 a = 1; a <= max_a; a++) {
                    i64 aa = (a * a) % MOD;
                    powA[a] = (aa * a) % MOD;
                }
            } else {
                for (i64 a = 1; a <= max_a; a++)
                    powA[a] = mod_pow(a, L);
            }

            for (i64 q = 0; q < q_full; q++) {
                i64 A = q + 1;
                i64 B = q + 2;
                i64 A_L = powA[A];
                i64 B_L = powA[B];
                i64 A_L1 = (i64)((i128)A_L * A % MOD);
                i128 term128 = (i128)q * A_L + (i128)(A * A % MOD) * B_L - (i128)B * A_L1;
                i64 term = (i64)(term128 % MOD);
                if (term < 0) term += MOD;
                total = (total + term) % MOD;
            }

            /* Last partial block for q = q_full */
            i64 q = q_full;
            m = (R - 1) - q_full * L;
            i64 A = q + 1;
            i64 B = q + 2;
            i64 A_L = powA[A];
            i64 term = (i64)((i128)q * A_L % MOD);
            if (m >= 1) {
                i64 A_L1 = (i64)((i128)A_L * A % MOD);
                i64 exp = L + 1 - m;
                i64 A_L1_m, B_m;

                if (exp == 1) A_L1_m = A % MOD;
                else if (exp == 2) A_L1_m = (A * A) % MOD;
                else if (exp == 3) A_L1_m = (i64)((i128)(A * A % MOD) * A % MOD);
                else A_L1_m = mod_pow(A, exp);

                if (m == 1) B_m = B % MOD;
                else if (m == 2) B_m = (B * B) % MOD;
                else if (m == 3) B_m = (i64)((i128)(B * B % MOD) * B % MOD);
                else B_m = mod_pow(B, m);

                i64 inner = (i64)(((i128)A_L1_m * B_m % MOD - A_L1) % MOD);
                if (inner < 0) inner += MOD;
                term = (i64)((i128)term + (i128)B * inner % MOD) % MOD;
            }
            total = (total + term) % MOD;
        }

        /* mu = 0 contribution for rem = R */
        i64 q = R / L;
        i64 r = R - q * L;
        i64 A = q + 1;
        i64 B = q + 2;
        i64 base;
        if (r == 0) {
            base = mod_pow(A, L);
        } else {
            base = (i64)((i128)mod_pow(A, L - r) * mod_pow(B, r) % MOD);
        }
        total = (total + base) % MOD;
    }

    free(powA);
    return total % MOD;
}
