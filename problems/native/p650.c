#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;
enum { LIMIT = 20000, MOD = 1000000007LL };

static i64 mod_pow(i64 a, i64 e) {
    i64 r = 1; a %= MOD;
    while (e > 0) {
        if (e & 1) r = (i64)((__int128)r * a % MOD);
        a = (i64)((__int128)a * a % MOD);
        e >>= 1;
    }
    return r;
}

long long pe650_answer(void) {
    int *spf = calloc(LIMIT + 1, sizeof(int));
    int *primes = calloc(LIMIT, sizeof(int));
    int pc = 0;
    for (int x = 2; x <= LIMIT; x++) {
        if (!spf[x]) { spf[x] = x; primes[pc++] = x; }
        for (int i = 0; i < pc; i++) {
            long y = (long)primes[i] * x;
            if (y > LIMIT) break;
            spf[y] = primes[i];
            if (primes[i] == spf[x]) break;
        }
    }
    int *pidx = calloc(LIMIT + 1, sizeof(int));
    for (int i = 0; i < pc; i++) pidx[primes[i]] = i;
    i64 *prime_power = calloc(pc, sizeof(i64));
    i64 *inv_fact = calloc(pc, sizeof(i64));
    i64 *inv_prime = calloc(pc, sizeof(i64));
    i64 *sigma_den = calloc(pc, sizeof(i64));
    for (int i = 0; i < pc; i++) {
        prime_power[i] = primes[i] % MOD;
        inv_fact[i] = 1;
        inv_prime[i] = mod_pow(primes[i], MOD - 2);
        sigma_den[i] = mod_pow(primes[i] - 1, MOD - 2);
    }
    int active = 0;
    i64 total = 1;
    for (int n = 2; n <= LIMIT; n++) {
        while (active < pc && primes[active] <= n) active++;
        for (int i = 0; i < active; i++)
            prime_power[i] = (i64)((__int128)prime_power[i] * inv_fact[i] % MOD);
        int x = n;
        while (x > 1) {
            int p = spf[x], e = 0;
            while (x % p == 0) { x /= p; e++; }
            int i = pidx[p];
            prime_power[i] = (i64)((__int128)prime_power[i] * mod_pow(p, (i64)(n - 1) * e) % MOD);
            inv_fact[i] = (i64)((__int128)inv_fact[i] * mod_pow(inv_prime[i], e) % MOD);
        }
        i64 divisor_sum = 1;
        for (int i = 0; i < active; i++) {
            i64 t = (prime_power[i] - 1) % MOD; if (t < 0) t += MOD;
            divisor_sum = (i64)((__int128)divisor_sum * t % MOD * sigma_den[i] % MOD);
        }
        total = (total + divisor_sum) % MOD;
    }
    free(spf); free(primes); free(pidx);
    free(prime_power); free(inv_fact); free(inv_prime); free(sigma_den);
    return total;
}
