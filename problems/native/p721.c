#include <math.h>
#include <stdint.h>
typedef long long i64;
enum { MOD = 999999937LL };
static void bin_exp(i64 x, i64 a, i64 n, i64 *a_sq, i64 *b_sq) {
    i64 A = x, B = 1;
    /* process bits of n after leading 1 */
    int msb = 63 - __builtin_clzll((unsigned long long)n);
    for (int i = msb - 1; i >= 0; i--) {
        i64 na = (A * A + a % MOD * (B % MOD) % MOD * (B % MOD)) % MOD;
        i64 nb = (2 * A % MOD * B) % MOD;
        A = na; B = nb;
        if ((n >> i) & 1) {
            na = (x % MOD * A + a % MOD * B) % MOD;
            nb = (A + x % MOD * B) % MOD;
            A = na; B = nb;
        }
    }
    *a_sq = A; *b_sq = B;
}
long long pe_solve(void) {
    i64 n = 5000000;
    i64 total = 0;
    for (i64 a = 1; a <= n; a++) {
        i64 ca = (i64)sqrt((double)a);
        while ((ca + 1) * (ca + 1) <= a) ca++;
        while (ca * ca > a) ca--;
        i64 asq, bsq;
        if (ca * ca == a) {
            bin_exp(ca, a, a * a, &asq, &bsq);
            total += 2 * asq;
        } else {
            bin_exp(ca + 1, a, a * a, &asq, &bsq);
            total += 2 * asq - 1;
        }
        total %= MOD;
    }
    return total;
}
