/* Project Euler 952: Order Modulo Factorial
 * Compute R(10^9+7, 10^7) mod 10^9+7.
 * Multiplicative order of prime p modulo n!.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef uint32_t u32;
typedef int64_t i64;
typedef uint64_t u64;
typedef __int128 i128;

#define P 1000000007LL
#define N 10000000

static u32 *spf;
static u32 *primes;
static u32 num_primes;
static u32 *max_exp;

static i64 mod_pow_u64(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    a %= mod;
    while (e > 0) {
        if (e & 1) r = (u64)((__int128)r * a % mod);
        a = (u64)((__int128)a * a % mod);
        e >>= 1;
    }
    return (i64)r;
}

static void linear_sieve_spf(u32 limit) {
    spf = (u32*)calloc(limit + 1, sizeof(u32));
    primes = (u32*)malloc(sizeof(u32) * 700000); /* pi(10^7) ~ 620k */

    num_primes = 0;
    for (u32 i = 2; i <= limit; i++) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes[num_primes++] = i;
        }
        for (u32 j = 0; j < num_primes; j++) {
            u32 p = primes[j];
            u32 ip = i * p;
            if (ip > limit) break;
            spf[ip] = p;
            if (p == spf[i]) break;
        }
    }
}

static u32 v_p_factorial(u32 n, u32 p) {
    u32 s = 0;
    while (n) {
        n /= p;
        s += n;
    }
    return s;
}

static u32 v2_val(u64 x) {
    u32 c = 0;
    while ((x & 1) == 0) {
        x >>= 1;
        c++;
    }
    return c;
}

static u32 order_mod_2_power_exponent(u64 a, u32 k) {
    if (k <= 1) return 0;
    if (k == 2) return ((a & 3) == 1) ? 0 : 1;

    if ((a & 3) == 1) {
        /* v2(a^{2^t}-1) = v2(a-1) + t */
        i64 t = (i64)k - (i64)v2_val(a - 1);
        return t > 0 ? (u32)t : 0;
    } else {
        /* a ≡ 3 mod 4: v2(a^{2^t}-1) = v2(a-1)+v2(a+1)+t-1 for t>=1 */
        i64 t = (i64)k - (i64)v2_val(a + 1);
        return t > 1 ? (u32)t : 1;
    }
}

/* Factorize x using SPF. Returns number of (prime, exponent) pairs. */
static u32 factorize(u32 x, u32 *fac_p, u32 *fac_e) {
    u32 cnt = 0;
    while (x > 1) {
        u32 p = spf[x];
        u32 e = 0;
        while (x % p == 0) {
            x /= p;
            e++;
        }
        fac_p[cnt] = p;
        fac_e[cnt] = e;
        cnt++;
    }
    return cnt;
}

static u64 multiplicative_order_mod_prime(u32 p, u32 q, u32 *max_exp_arr) {
    /* Compute ord_q(p) and update max_exp for prime factors of ord */
    u32 m = q - 1;
    u32 fac_p[32], fac_e[32];
    u32 nfac = factorize(m, fac_p, fac_e);

    u64 r = m;
    u64 base = p % q;

    for (u32 i = 0; i < nfac; i++) {
        u32 f = fac_p[i];
        u32 e_rem = fac_e[i];
        while (e_rem) {
            u64 cand = r / f;
            if (mod_pow_u64(base, cand, q) == 1) {
                r = cand;
                e_rem--;
            } else {
                break;
            }
        }
        if (e_rem > max_exp_arr[f]) {
            max_exp_arr[f] = e_rem;
        }
    }
    return r;
}

static u32 q_adic_valuation(u64 p, u64 r, u32 q, u32 limit) {
    if (limit <= 1) return 1;
    u32 s = 1;
    u64 q_pow = q;
    while (s < limit) {
        u64 q_pow_next = q_pow * q;
        /* Check p^r mod q^(s+1) == 1 */
        if (mod_pow_u64(p, r, q_pow_next) != 1) break;
        q_pow = q_pow_next;
        s++;
    }
    return s;
}

long long p952_native(void) {
    linear_sieve_spf(N);
    max_exp = (u32*)calloc(N + 1, sizeof(u32));

    /* Handle q=2 separately */
    u32 a2 = v_p_factorial(N, 2);
    u32 t2 = order_mod_2_power_exponent((u64)P, a2);
    if (t2 > max_exp[2]) max_exp[2] = t2;

    /* Odd primes */
    for (u32 qi = 0; qi < num_primes; qi++) {
        u32 q = primes[qi];
        if (q == 2) continue;

        u32 a = v_p_factorial(N, q);
        u64 r0 = multiplicative_order_mod_prime((u32)P, q, max_exp);
        u32 s = q_adic_valuation((u64)P, r0, q, a);
        i64 extra = (i64)a - (i64)s;
        if (extra > 0 && (u32)extra > max_exp[q]) {
            max_exp[q] = (u32)extra;
        }
    }

    /* Reconstruct order modulo P */
    u64 res = 1 % (u64)P;
    for (u32 qi = 0; qi < num_primes; qi++) {
        u32 q = primes[qi];
        u32 e = max_exp[q];
        if (e) {
            res = (u64)((__int128)res * mod_pow_u64(q, e, (u64)P) % (u64)P);
        }
    }

    free(spf);
    free(primes);
    free(max_exp);

    return (long long)res;
}
