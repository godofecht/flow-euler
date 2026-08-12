#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef int64_t i64;

/* ---------- Sieve ---------- */
static char *is_comp;       /* 0 = prime, 1 = composite, for 0..N */
static int *primes_arr;
static int num_primes;

static void sieve_upto(int n) {
    is_comp = (char *)calloc(n + 1, 1);
    is_comp[0] = is_comp[1] = 1;
    int limit = (int)sqrtl((long double)n);
    for (int i = 2; i <= limit; i++) {
        if (!is_comp[i]) {
            for (int j = i * i; j <= n; j += i)
                is_comp[j] = 1;
        }
    }
    /* count primes */
    int cnt = 0;
    for (int i = 2; i <= n; i++) if (!is_comp[i]) cnt++;
    primes_arr = (int *)malloc(cnt * sizeof(int));
    num_primes = 0;
    for (int i = 2; i <= n; i++) if (!is_comp[i]) primes_arr[num_primes++] = i;
}

/* ---------- Modular exponentiation ---------- */
static i64 mod_pow_i64(i64 a, i64 e, i64 mod) {
    i64 r = 1 % mod; a %= mod;
    if (a < 0) a += mod;
    while (e > 0) {
        if (e & 1) r = r * a % mod;
        a = a * a % mod;
        e >>= 1;
    }
    return r;
}

/* ---------- Cycle detection: does orbit of 1 under x->1+x^k hit 0 mod q? ---------- */
static int hits_zero_mod_prime(i64 q, i64 k) {
    if (q == 2) return 1;
    if ((q - 1) % k != 0) return 1;  /* gcd(k,q-1)=1 => permutation => must hit 0 */

    /* f(x) = (1 + x^k) % q */
    /* Brent's cycle detection */
    i64 tortoise = 1;
    i64 hare = (1 + mod_pow_i64(1, k, q)) % q;  /* f(1) = 1 + 1 = 2, but use general */
    if (hare == 0) return 1;

    i64 power = 1, lam = 1;
    while (tortoise != hare) {
        if (hare == 0) return 1;
        if (power == lam) {
            tortoise = hare;
            power <<= 1;
            lam = 0;
        }
        hare = (1 + mod_pow_i64(hare, k, q)) % q;
        lam++;
    }
    return hare == 0;
}

/* For k=2, f(x) = x*x + 1 mod q.  Specialized for speed. */
static int hits_zero_mod_prime_k2(i64 q) {
    if (q == 2) return 1;
    /* (q-1) % 2 == 0 always for odd q, so we always need cycle detection */
    i64 tortoise = 1;
    i64 hare = (1 + 1) % q;  /* 1*1 + 1 = 2 */
    if (hare == 0) return 1;

    i64 power = 1, lam = 1;
    while (tortoise != hare) {
        if (hare == 0) return 1;
        if (power == lam) {
            tortoise = hare;
            power <<= 1;
            lam = 0;
        }
        hare = (hare * hare + 1) % q;
        lam++;
    }
    return hare == 0;
}

/* ---------- Distinct prime factors of n (n <= 10^7) ---------- */
static int prime_factors_unique(int n, int *out) {
    int x = n, cnt = 0;
    for (int i = 0; i < num_primes; i++) {
        int p = primes_arr[i];
        if ((i64)p * p > x) break;
        if (x % p == 0) {
            out[cnt++] = p;
            while (x % p == 0) x /= p;
        }
    }
    if (x > 1) out[cnt++] = x;
    return cnt;
}

/* ---------- Find good primes ---------- */
static int *good_primes;
static int num_good;

static void find_good_primes(int limit) {
    good_primes = (int *)malloc(num_primes * sizeof(int));
    num_good = 0;

    int factors_buf[32];

    for (int pi = 0; pi < num_primes; pi++) {
        int q = primes_arr[pi];
        if (q == 2) {
            good_primes[num_good++] = q;
            continue;
        }
        /* odd q: must be 1 mod 4 */
        if ((q & 3) != 1) continue;
        /* check exponent 2 */
        if (!hits_zero_mod_prime_k2(q)) continue;

        /* check each odd prime factor of q-1 */
        int nf = prime_factors_unique(q - 1, factors_buf);
        int ok = 1;
        for (int fi = 0; fi < nf; fi++) {
            int p = factors_buf[fi];
            if (p == 2) continue;
            if (!hits_zero_mod_prime(q, p)) {
                ok = 0;
                break;
            }
        }
        if (ok) good_primes[num_good++] = q;
    }
}

/* ---------- Sum of squarefree products <= limit ---------- */
static i64 sum_squarefree_products_leq(int limit, int *plist, int plen) {
    /* dynamic array of products */
    i64 cap = 1 << 20;
    int *prods = (int *)malloc(cap * sizeof(int));
    int len = 1;
    prods[0] = 1;

    for (int pi = 0; pi < plen; pi++) {
        int p = plist[pi];
        int base_len = len;
        for (int i = 0; i < base_len; i++) {
            i64 v = (i64)prods[i] * p;
            if (v <= limit) {
                if (len >= cap) {
                    cap <<= 1;
                    prods = (int *)realloc(prods, cap * sizeof(int));
                }
                prods[len++] = (int)v;
            }
        }
    }

    i64 total = 0;
    for (int i = 0; i < len; i++) total += prods[i];
    free(prods);
    return total;
}

long long p927_native(void) {
    int limit = 10000000;
    sieve_upto(limit);
    find_good_primes(limit);
    i64 result = sum_squarefree_products_leq(limit, good_primes, num_good);
    free(is_comp);
    free(primes_arr);
    free(good_primes);
    return (long long)result;
}
