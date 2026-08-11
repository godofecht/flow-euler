#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef int64_t i64;

#define MOD 999999001LL
#define ROOT 17LL
#define N 1000

static i64 mod_pow(i64 a, i64 e, i64 mod) {
    i64 r = 1 % mod; a %= mod; if (a < 0) a += mod;
    while (e > 0) {
        if (e & 1) r = (i64)((__int128)r * a % mod);
        a = (i64)((__int128)a * a % mod);
        e >>= 1;
    }
    return r;
}

static i64 mod_inv(i64 a, i64 mod) { return mod_pow(a, mod - 2, mod); }

static int sieve_primes(int N_lim, i64 *primes) {
    char *is_prime = (char *)malloc((N_lim + 1) * sizeof(char));
    memset(is_prime, 1, (N_lim + 1) * sizeof(char));
    is_prime[0] = is_prime[1] = 0;
    int count = 0;
    for (int i = 2; i <= N_lim; i++) {
        if (is_prime[i]) {
            primes[count++] = i;
            for (long long j = (long long)i * i; j <= N_lim; j += i)
                is_prime[j] = 0;
        }
    }
    free(is_prime);
    return count;
}

/* For each prime p, compute E = sum_{k=1}^{N} (N+1-k) * v_p(k!). */
static void compute_exponents(int N_lim, const i64 *primes, int np, i64 *exps) {
    int *vp_in_k = (int *)malloc((N_lim + 1) * sizeof(int));
    i64 *fact_vp = (i64 *)malloc((N_lim + 1) * sizeof(i64));
    for (int pi = 0; pi < np; pi++) {
        i64 p = primes[pi];
        for (int k = 0; k <= N_lim; k++) vp_in_k[k] = 0;
        for (int k = 1; k <= N_lim; k++) {
            int x = k, c = 0;
            while (x % p == 0) { x /= (int)p; c++; }
            vp_in_k[k] = c;
        }
        fact_vp[0] = 0;
        for (int k = 1; k <= N_lim; k++)
            fact_vp[k] = fact_vp[k - 1] + vp_in_k[k];
        i64 E = 0;
        for (int k = 1; k <= N_lim; k++)
            E += (i64)(N_lim + 1 - k) * fact_vp[k];
        exps[pi] = E;
    }
    free(vp_in_k);
    free(fact_vp);
}

static i64 D_star_mod(int N_lim) {
    i64 primes[256];
    int np = sieve_primes(N_lim, primes);
    i64 *exps = (i64 *)malloc(np * sizeof(i64));
    compute_exponents(N_lim, primes, np, exps);

    i64 m = N_lim;
    i64 omega = mod_pow(ROOT, (MOD - 1) / m, MOD);
    i64 inv_m = mod_inv(m, MOD);
    i64 ans = 0;
    i64 w = 1;
    for (i64 t = 0; t < m; t++) {
        i64 prod = 1;
        for (int pi = 0; pi < np; pi++) {
            i64 p = primes[pi];
            i64 E = exps[pi];
            i64 r = (i64)((__int128)p * w % MOD);
            i64 term;
            if (r == 1) {
                term = (E + 1) % MOD;
            } else {
                i64 num = (mod_pow(r, E + 1, MOD) - 1 + MOD) % MOD;
                i64 den = (r - 1 + MOD) % MOD;
                i64 den_inv = mod_inv(den, MOD);
                term = (i64)((__int128)num * den_inv % MOD);
            }
            prod = (i64)((__int128)prod * term % MOD);
        }
        ans = (ans + prod) % MOD;
        w = (i64)((__int128)w * omega % MOD);
    }
    i64 result = (i64)((__int128)ans * inv_m % MOD);
    free(exps);
    return result;
}

long long p956_native(void) {
    return (long long)D_star_mod(N);
}
