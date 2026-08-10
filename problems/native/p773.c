#include <stdint.h>
#include <stdlib.h>
typedef long long i64;
enum { MOD = 1000000007LL };
static i64 modpow(i64 b, i64 e) {
    i64 r = 1; b %= MOD;
    while (e > 0) { if (e & 1) r = r * b % MOD; b = b * b % MOD; e >>= 1; }
    return r;
}
static int is_prime(int n) {
    if (n < 2) return 0;
    for (int i = 2; i * i <= n; i++) if (n % i == 0) return 0;
    return 1;
}
long long pe_solve(void) {
    int k = 97;
    int *primes = malloc((size_t)k * sizeof(int));
    int cnt = 0, x = 7;
    while (cnt < k) {
        if (x % 10 == 7 && is_prime(x)) primes[cnt++] = x;
        x += 2;
    }
    i64 M_mod = 1, phi_mod = 1;
    for (int i = 0; i < k; i++) {
        M_mod = M_mod * primes[i] % MOD;
        phi_mod = phi_mod * (primes[i] - 1) % MOD;
    }
    int q_table[4] = {7, 1, 3, 9};
    i64 A = 0, c = 1;
    for (int s = 0; s <= k; s++) {
        i64 term = c * q_table[s & 3] % MOD;
        if (s & 1) A = (A - term) % MOD;
        else A = (A + term) % MOD;
        if (s < k) {
            c = c * (k - s) % MOD;
            c = c * modpow(s + 1, MOD - 2) % MOD;
        }
    }
    if (A < 0) A += MOD;
    free(primes);
    return M_mod * ((A + 5 * phi_mod) % MOD) % MOD;
}
