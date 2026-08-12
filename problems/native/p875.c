// Project Euler 875: Quadruple Congruence.
// Q(N) = sum_{n=1..N} q(n) mod 1001961001, where q(n) counts solutions to
// a1^2+a2^2+a3^2+a4^2 = b1^2+b2^2+b3^2+b4^2 (mod n).
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned long long u64;
typedef uint32_t u32;
typedef __int128 i128;

#define MOD 1001961001ULL
#define TARGET_N 12345678

static u64 mulmod(u64 a, u64 b, u64 m) {
    return (u64)((i128)a * b % m);
}

static u64 powmod(u64 a, u64 e, u64 m) {
    u64 r = 1 % m;
    a %= m;
    while (e > 0) {
        if (e & 1) r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}

static u64 q_prime_power_mod(u64 p, int k, u64 mod) {
    if (p == 2) {
        if (k == 1) return 128 % mod;
        u64 r = 8 % mod;
        u64 g = 0, cur = 1;
        for (int i = 0; i < k - 1; i++) {
            g += cur;
            if (g >= mod) g -= mod;
            cur = mulmod(cur, r, mod);
        }
        u64 term1 = powmod(2, 7ULL * k, mod);
        u64 term2 = mulmod(powmod(2, 4ULL * k + 3, mod), g, mod);
        return (term1 + term2) % mod;
    }
    // Odd prime
    if (k == 1) {
        u64 p2 = mulmod(p, p, mod);
        u64 p3 = mulmod(p2, p, mod);
        u64 p4 = mulmod(p2, p2, mod);
        u64 p7 = mulmod(p4, p3, mod);
        return (p7 + p4 + mod - p3) % mod;
    }
    u64 r = powmod(p, 3, mod);
    u64 g = 0, cur = 1;
    for (int i = 0; i < k; i++) {
        g += cur;
        if (g >= mod) g -= mod;
        cur = mulmod(cur, r, mod);
    }
    u64 term1 = powmod(p, 7ULL * k, mod);
    u64 termp = powmod(p, 4ULL * k - 1, mod);
    u64 term2 = mulmod(mulmod((p - 1) % mod, termp, mod), g, mod);
    return (term1 + term2) % mod;
}

// Linear sieve for smallest prime factor.
static void smallest_prime_factors(u64 n, u32 *spf, u32 *primes, int *prime_count) {
    memset(spf, 0, (n + 1) * sizeof(u32));
    spf[1] = 1;
    int pc = 0;
    for (u64 i = 2; i <= n; i++) {
        if (spf[i] == 0) {
            spf[i] = (u32)i;
            primes[pc++] = (u32)i;
        }
        u64 si = spf[i];
        for (int j = 0; j < pc; j++) {
            u64 p = primes[j];
            u64 ip = i * p;
            if (ip > n) break;
            spf[ip] = (u32)p;
            if (p == si) break;
        }
    }
    *prime_count = pc;
}

static u64 compute_Q(u64 N, u64 mod) {
    u32 *spf = (u32 *)malloc((N + 1) * sizeof(u32));
    u32 *primes = (u32 *)malloc((N + 1) * sizeof(u32));
    int pc;
    smallest_prime_factors(N, spf, primes, &pc);

    u32 *f = (u32 *)calloc(N + 1, sizeof(u32));
    f[1] = 1;
    u64 total = 1 % mod;

    for (u64 n = 2; n <= N; n++) {
        u64 p = spf[n];
        u64 m = n / p;
        u64 val;

        if (spf[m] != p) {
            u64 qp;
            if (m == 1) {
                qp = q_prime_power_mod(p, 1, mod);
            } else {
                qp = f[p];
            }
            val = mulmod(f[m], qp, mod);
        } else {
            u64 ppow = p * p;
            int k = 2;
            u64 mm = m / p;
            while (mm > 1 && spf[mm] == p) {
                ppow *= p;
                k++;
                mm /= p;
            }
            u64 rest = mm;
            u64 qp;
            if (rest == 1) {
                qp = q_prime_power_mod(p, k, mod);
            } else {
                qp = f[ppow];
            }
            val = mulmod(f[rest], qp, mod);
        }

        f[n] = (u32)val;
        total += val;
        if (total >= mod) total -= mod;
    }

    free(spf);
    free(primes);
    free(f);
    return total;
}

long long p875_native(void) {
    return (long long)compute_Q(TARGET_N, MOD);
}
