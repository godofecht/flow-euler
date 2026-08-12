// Project Euler 907: Stacking Cups
// S(n) satisfies an order-8 linear recurrence for n >= 10:
// S(n) = 2S(n-1) - 3S(n-2) + 5S(n-3) - 4S(n-4) + 4S(n-5) - 3S(n-6) + S(n-7) - S(n-8)
// Compute S(10^7) mod 1e9+7 via matrix exponentiation.
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define MOD 1000000007LL

// Exact small values S(1)..S(9)
static const long long small[10] = {0, 2, 2, 6, 12, 16, 22, 36, 58, 82};

// Recurrence coefficients (for S(n) in terms of S(n-1)..S(n-8))
static const long long coeff[8] = {2, -3, 5, -4, 4, -3, 1, -1};

static void mat_mul(long long A[8][8], long long B[8][8], long long C[8][8]) {
    long long tmp[8][8];
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            long long s = 0;
            for (int k = 0; k < 8; k++) {
                s = (s + A[i][k] * B[k][j]) % MOD;
            }
            C[i][j] = (s % MOD + MOD) % MOD;
        }
    }
}

static void mat_vec_mul(long long A[8][8], long long v[8], long long out[8]) {
    for (int i = 0; i < 8; i++) {
        long long s = 0;
        for (int j = 0; j < 8; j++) {
            s = (s + A[i][j] * v[j]) % MOD;
        }
        out[i] = (s % MOD + MOD) % MOD;
    }
}

long long p907_native(void) {
    long long n = 10000000LL;
    if (n <= 9) return small[n] % MOD;

    // Build companion matrix
    long long M[8][8];
    memset(M, 0, sizeof(M));
    for (int j = 0; j < 8; j++) {
        M[0][j] = ((coeff[j] % MOD) + MOD) % MOD;
    }
    for (int i = 1; i < 8; i++) {
        M[i][i - 1] = 1;
    }

    // Base state at t=9: [S(9), S(8), ..., S(2)]
    long long base[8] = {
        small[9] % MOD, small[8] % MOD, small[7] % MOD, small[6] % MOD,
        small[5] % MOD, small[4] % MOD, small[3] % MOD, small[2] % MOD
    };

    // Compute M^(n-9) * base
    long long e = n - 9;
    long long result[8][8];
    // identity
    memset(result, 0, sizeof(result));
    for (int i = 0; i < 8; i++) result[i][i] = 1;

    long long base_pow[8][8];
    memcpy(base_pow, M, sizeof(M));

    while (e > 0) {
        if (e & 1) {
            long long tmp[8][8];
            mat_mul(result, base_pow, tmp);
            memcpy(result, tmp, sizeof(tmp));
        }
        e >>= 1;
        if (e > 0) {
            long long tmp[8][8];
            mat_mul(base_pow, base_pow, tmp);
            memcpy(base_pow, tmp, sizeof(tmp));
        }
    }

    long long out[8];
    mat_vec_mul(result, base, out);
    return out[0];
}
