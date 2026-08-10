#include <stdint.h>
typedef long long i64;
enum { MOD = 1000000007LL };
static i64 modpow(i64 b, i64 e) {
    i64 r = 1;
    b %= MOD;
    while (e > 0) {
        if (e & 1) r = r * b % MOD;
        b = b * b % MOD;
        e >>= 1;
    }
    return r;
}
static i64 sum_pows_8(i64 k) {
    if (k <= 0) return 0;
    return (modpow(8, k) - 1 + MOD) % MOD * modpow(7, MOD - 2) % MOD;
}
static i64 sum_pows_4_1_to_m(i64 m) {
    if (m <= 0) return 0;
    return (modpow(4, m + 1) - 4 + MOD) % MOD * modpow(3, MOD - 2) % MOD;
}
static i64 sum_T_upto(i64 K) {
    if (K <= 0) return 0;
    i64 T = 0, s = 0, pow2 = 1, pow4 = 1;
    for (i64 i = 0; i < K; i++) {
        i64 add = (pow2 + 2) % MOD * pow4 % MOD;
        T = (T + T + add) % MOD;
        s += T; if (s >= MOD) s -= MOD;
        pow2 = (pow2 + pow2) % MOD;
        pow4 = pow4 * 4 % MOD;
    }
    return s;
}
static i64 S(i64 N) {
    if (N % 2 == 0) {
        i64 m = N / 2;
        i64 A = (sum_pows_8(m) + sum_T_upto(m - 1)) % MOD;
        i64 B = (sum_pows_4_1_to_m(m) - m % MOD + MOD) % MOD;
        return (A + modpow(4, m) + B) % MOD;
    } else {
        i64 m = (N - 1) / 2;
        i64 A = (sum_pows_8(m + 1) + sum_T_upto(m)) % MOD;
        i64 B = (sum_pows_4_1_to_m(m) - m % MOD + MOD) % MOD;
        return (A + B) % MOD;
    }
}
long long pe_solve(void) { return S(12345678); }
