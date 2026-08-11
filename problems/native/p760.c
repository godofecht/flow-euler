#include <stdint.h>
#include <string.h>
typedef long long i64;
enum { MOD = 1000000007LL };

/*
 * Project Euler 760: Sum over Bitwise Operators.
 *
 * g(m,n) = (m xor n) + (m or n) + (m and n) = 2*(m|n).
 * G(N) = sum_{n=0..N} sum_{k=0..n} g(k, n-k).
 * Answer: G(10^18) mod 1e9+7.
 *
 * Sum of (a|b) over all pairs (a,b) with a+b<=N equals
 *   sum_i 2^i * (#pairs where OR-bit i is 1).
 * OR-bit i is 0 iff a_i=0 and b_i=0, counted via a digit-DP
 * over the binary sum a+b enforcing S<=N with both bits forced to 0.
 */

static i64 count_pairs_sum_leq(i64 N, int fixed_bit) {
    if (N < 0) return 0;
    int bits = 1;
    for (i64 t = N; t > 1; t >>= 1) bits++;

    /* dp[carry_next][less] */
    i64 dp[2][2], ndp[2][2];
    memset(dp, 0, sizeof(dp));
    dp[0][0] = 1;

    for (int pos = bits - 1; pos >= 0; pos--) {
        int nbit = (int)((N >> pos) & 1);
        memset(ndp, 0, sizeof(ndp));

        for (int carry_next = 0; carry_next < 2; carry_next++) {
            for (int less = 0; less < 2; less++) {
                i64 ways = dp[carry_next][less];
                if (ways == 0) continue;

                for (int carry_cur = 0; carry_cur < 2; carry_cur++) {
                    for (int a_bit = 0; a_bit < 2; a_bit++) {
                        for (int b_bit = 0; b_bit < 2; b_bit++) {
                            if (fixed_bit >= 0 && pos == fixed_bit && (a_bit | b_bit))
                                continue;

                            int total = a_bit + b_bit + carry_cur;
                            if ((total >> 1) != carry_next) continue;

                            int s_bit = total & 1;
                            if (less == 0 && s_bit > nbit) continue;

                            int new_less = less | (s_bit < nbit);
                            ndp[carry_cur][new_less] =
                                (ndp[carry_cur][new_less] + ways) % MOD;
                        }
                    }
                }
            }
        }
        memcpy(dp, ndp, sizeof(dp));
    }
    /* carry into the bit below LSB must be 0 */
    return (dp[0][0] + dp[0][1]) % MOD;
}

long long p760_native(void) {
    i64 N = 1000000000000000000LL; /* 10^18 */

    i64 inv2 = (MOD + 1) / 2;
    i64 total_pairs = ((N + 1) % MOD) * ((N + 2) % MOD) % MOD;
    total_pairs = total_pairs * inv2 % MOD;

    int bits = 1;
    for (i64 t = N; t > 1; t >>= 1) bits++;

    i64 pow2 = 1, sum_or = 0;
    for (int i = 0; i < bits; i++) {
        i64 both_zero = count_pairs_sum_leq(N, i);
        i64 bit_is_one = (total_pairs - both_zero) % MOD;
        if (bit_is_one < 0) bit_is_one += MOD;
        sum_or = (sum_or + pow2 * bit_is_one) % MOD;
        pow2 = (pow2 * 2) % MOD;
    }

    return (2 * sum_or) % MOD;
}
