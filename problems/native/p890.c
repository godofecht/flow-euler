// Project Euler 890: Binary Partitions.
// p(n) = number of partitions of n into powers of 2.
// p(2m) = p(2m+1) = S(m) = [x^m] A(x) where A(x) = prod_{k>=0} (1+x^{2^k})^{k+2}.
// Carry-DP in base 2 with convolution via GMP big-integer multiplication.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

typedef unsigned long long u64;
typedef __int128 i128;

#define MOD 1000000007ULL

static u64 mulmod(u64 a, u64 b, u64 m) {
    return (u64)((i128)a * b % m);
}

static u64 powmod(u64 a, u64 e, u64 m) {
    u64 r = 1 % m; a %= m;
    while (e > 0) {
        if (e & 1) r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}

#define BASE_BITS 80

// Convolve a and b (integer convolution), then decimate: take indices bit, bit+2, ...
// Results reduced mod MOD. Uses GMP for fast multiplication.
static void convolve_and_decimate(const u64 *a, int la, const u64 *b, int lb,
                                   int bit, u64 *res, int *res_len) {
    int out_len = la + lb - 1;

    mpz_t ia, ib, prod;
    mpz_init(ia);
    mpz_init(ib);
    mpz_init(prod);

    for (int i = la - 1; i >= 0; i--) {
        mpz_mul_2exp(ia, ia, BASE_BITS);
        mpz_add_ui(ia, ia, (unsigned long)a[i]);
    }
    for (int i = lb - 1; i >= 0; i--) {
        mpz_mul_2exp(ib, ib, BASE_BITS);
        mpz_add_ui(ib, ib, (unsigned long)b[i]);
    }

    mpz_mul(prod, ia, ib);

    *res_len = (out_len - bit + 1) / 2;

    mpz_t temp, digit, mod_val;
    mpz_init(temp);
    mpz_init(digit);
    mpz_init_set_ui(mod_val, MOD);

    // Shift right by bit * BASE_BITS to align to first needed digit
    if (bit > 0)
        mpz_fdiv_q_2exp(temp, prod, (mp_bitcnt_t)(bit * BASE_BITS));
    else
        mpz_set(temp, prod);

    for (int i = 0; i < *res_len; i++) {
        mpz_fdiv_r_2exp(digit, temp, BASE_BITS);
        mpz_mod(digit, digit, mod_val);
        res[i] = mpz_get_ui(digit);
        mpz_fdiv_q_2exp(temp, temp, (mp_bitcnt_t)(2 * BASE_BITS));
    }

    mpz_clear(ia);
    mpz_clear(ib);
    mpz_clear(prod);
    mpz_clear(temp);
    mpz_clear(digit);
    mpz_clear(mod_val);
}

static u64 fact[2500], invfact[2500];

static void prepare_factorials(int nmax) {
    fact[0] = 1;
    for (int i = 1; i <= nmax; i++)
        fact[i] = mulmod(fact[i - 1], (u64)i, MOD);
    invfact[nmax] = powmod(fact[nmax], MOD - 2, MOD);
    for (int i = nmax; i >= 1; i--)
        invfact[i - 1] = mulmod(invfact[i], (u64)i, MOD);
}

static u64 coefficient_A(const mpz_t m_val) {
    if (mpz_sgn(m_val) == 0) return 1;

    int L = (int)mpz_sizeinbase(m_val, 2);  // bit length of m
    int max_m = L + 2;
    prepare_factorials(max_m);

    int cap = L + 10;
    u64 *dp = (u64 *)malloc(cap * sizeof(u64));
    u64 *new_dp = (u64 *)malloc(cap * sizeof(u64));
    int dp_len = 1;
    dp[0] = 1;

    for (int k = 0; k < L; k++) {
        int bit = (int)mpz_tstbit(m_val, k);
        int top = k + 2;

        // Build binomial row C(top, j) mod MOD
        u64 row[2500];
        for (int j = 0; j <= top; j++)
            row[j] = mulmod(fact[top], mulmod(invfact[j], invfact[top - j], MOD), MOD);

        int new_len;
        convolve_and_decimate(dp, dp_len, row, top + 1, bit, new_dp, &new_len);

        // Swap dp and new_dp
        memcpy(dp, new_dp, new_len * sizeof(u64));
        dp_len = new_len;
    }

    u64 result = dp[0] % MOD;
    free(dp);
    free(new_dp);
    return result;
}

long long p890_native(void) {
    mpz_t seven, n_val, m_val;
    mpz_init(seven);
    mpz_init(n_val);
    mpz_init(m_val);

    // n = 7^777
    mpz_set_ui(seven, 7);
    mpz_pow_ui(n_val, seven, 777);

    // m = n // 2
    mpz_fdiv_q_2exp(m_val, n_val, 1);

    u64 ans = coefficient_A(m_val);

    mpz_clear(seven);
    mpz_clear(n_val);
    mpz_clear(m_val);

    return (long long)ans;
}
