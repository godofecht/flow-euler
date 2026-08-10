#include <stdint.h>
#include <stdlib.h>
typedef long long i64;
enum { MOD = 1000000007LL, BLOCK = 200000 };
static i64 modpow(i64 b, i64 e) {
    i64 r = 1; b %= MOD;
    while (e > 0) { if (e & 1) r = r * b % MOD; b = b * b % MOD; e >>= 1; }
    return r;
}
static void inv_range(i64 start, i64 end, i64 *out) {
    i64 L = end - start + 1;
    i64 prod = 1, x = start;
    for (i64 i = 0; i < L; i++) { prod = prod * x % MOD; out[i] = prod; x++; }
    i64 inv_prod = modpow(prod, MOD - 2), suffix = inv_prod;
    x = end;
    for (i64 i = L - 1; i >= 0; i--) {
        i64 prev = i ? out[i - 1] : 1;
        out[i] = suffix * prev % MOD;
        suffix = suffix * x % MOD;
        x--;
    }
}
long long pe_solve(void) {
    i64 k = 100000000LL;
    i64 n = 10000000000000000LL;
    i64 m = n / k;
    i64 base = modpow(2, n); /* 2^n mod MOD via binary exp on n */
    /* modpow only takes i64 e - good for n=1e16 */
    i64 pow2_2m = modpow(2, 2 * m);
    i64 r = modpow(pow2_2m, MOD - 2);
    i64 tmax = k / 2;
    i64 term = 1, acc = 1, a = k;
    i64 *invs = malloc((size_t)BLOCK * sizeof(i64));
    i64 start = 1;
    while (start <= tmax) {
        i64 end = start + BLOCK - 1;
        if (end > tmax) end = tmax;
        inv_range(start, end, invs);
        i64 L = end - start + 1;
        for (i64 i = 0; i < L; i++) {
            i64 inv = invs[i];
            i64 invsq_r = inv * inv % MOD * r % MOD;
            term = term * a % MOD * (a - 1) % MOD * invsq_r % MOD;
            acc += term;
            a -= 2;
        }
        start = end + 1;
    }
    free(invs);
    return base * (acc % MOD) % MOD;
}
