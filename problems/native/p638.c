/* Project Euler 638 - Weighted Lattice Paths
 * sum_{k=1..7} [2*10^k+k choose 10^k+k]_q (mod 1e9+7)
 * Gaussian binomial coefficient via product formula.
 */
#include <stdint.h>

static int64_t mod_pow(int64_t base, int64_t exp, int64_t mod) {
    int64_t r = 1 % mod;
    base %= mod;
    if (base < 0) base += mod;
    while (exp > 0) {
        if (exp & 1) r = (__int128)r * base % mod;
        base = (__int128)base * base % mod;
        exp >>= 1;
    }
    return r;
}

static int64_t gaussian_binom(int64_t m, int64_t n, int64_t q, int64_t mod) {
    if (n < 0 || n > m) return 0;
    if (n == 0 || n == m) return 1;
    if (q % mod == 1) {
        /* regular binomial */
        if (n > m - n) n = m - n;
        if (n == 0) return 1;
        int64_t num = 1, den = 1;
        int64_t start = m - n;
        for (int64_t i = 1; i <= n; i++) {
            num = (__int128)num * ((start + i) % mod) % mod;
            den = (__int128)den * (i % mod) % mod;
        }
        return (__int128)num * mod_pow(den, mod - 2, mod) % mod;
    }
    n = (n < m - n) ? n : (m - n);
    if (n == 0) return 1;
    int64_t b = m - n;
    int64_t num_prod = 1, den_prod = 1;
    int64_t qmod = q % mod;
    if (qmod < 0) qmod += mod;
    int64_t pow_den = qmod;
    int64_t pow_num = mod_pow(qmod, b + 1, mod);
    for (int64_t i = 0; i < n; i++) {
        num_prod = (__int128)num_prod * (mod + 1 - pow_num) % mod;
        den_prod = (__int128)den_prod * (mod + 1 - pow_den) % mod;
        pow_den = (__int128)pow_den * qmod % mod;
        pow_num = (__int128)pow_num * qmod % mod;
    }
    return (__int128)num_prod * mod_pow(den_prod, mod - 2, mod) % mod;
}

long long p638_native(void) {
    int64_t mod = 1000000007LL;
    int64_t ans = 0;
    for (int k = 1; k <= 7; k++) {
        int64_t n = 1;
        for (int j = 0; j < k; j++) n *= 10;
        n += k;
        ans = (ans + gaussian_binom(2 * n, n, k, mod)) % mod;
    }
    return (long long)ans;
}
