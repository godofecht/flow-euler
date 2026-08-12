/* Project Euler 937 - Equiproduct Partition.
   Computes G(10^8) mod 1_000_000_007.

   Approach: k! is "equiproduct" iff the parity of popcount of cumulative
   exponents of relevant primes (p%8 in {5,7}, inert in Z[sqrt(-2)]) plus
   the p=2 exponent parity is 0.

   We decompose the computation:
   1. Sieve all odd primes up to n.
   2. For each relevant prime p (and p=2), iterate over multiples k,
      track the running exponent E_p, and flip a delta bit at k whenever
      popcount(E_p) parity changes.
   3. Final pass: compute running global parity from delta bits, and sum
      k! mod MOD when the parity is 0.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define MOD 1000000007ULL

static uint64_t n_val;
static uint8_t *is_prime;      /* odd-prime sieve, byte per odd number */
static uint64_t *delta;        /* bit array: parity flip at position k */

static inline void delta_flip(uint64_t k) {
    delta[k >> 6] ^= (1ULL << (k & 63));
}

static inline int delta_get(uint64_t k) {
    return (delta[k >> 6] >> (k & 63)) & 1;
}

static void build_sieve(uint64_t n) {
    uint64_t size = n / 2 + 1;
    is_prime = malloc(size);
    memset(is_prime, 1, size);
    is_prime[0] = 0;  /* 1 is not prime */
    uint64_t limit = (uint64_t)sqrt((double)n);
    for (uint64_t i = 1; 2 * i + 1 <= limit; i++) {
        if (is_prime[i]) {
            uint64_t p = 2 * i + 1;
            for (uint64_t j = p * p; j <= n; j += 2 * p) {
                is_prime[j / 2] = 0;
            }
        }
    }
}

static void process_prime_2(uint64_t n) {
    uint64_t E = 0;
    int par = 0;
    for (uint64_t k = 2; k <= n; k += 2) {
        int e = __builtin_ctzll(k);
        E += e;
        int new_par = __builtin_popcountll(E) & 1;
        if (new_par != par) {
            delta_flip(k);
            par = new_par;
        }
    }
}

static void process_relevant_prime(uint64_t p, uint64_t n) {
    uint64_t E = 0;
    int par = 0;
    for (uint64_t m = 1; m * p <= n; m++) {
        uint64_t k = m * p;
        int e = 1;
        uint64_t r = m;
        while (r % p == 0) {
            e++;
            r /= p;
        }
        E += e;
        int new_par = __builtin_popcountll(E) & 1;
        if (new_par != par) {
            delta_flip(k);
            par = new_par;
        }
    }
}

static uint64_t compute_G(uint64_t n) {
    n_val = n;
    build_sieve(n);

    uint64_t delta_words = (n + 64) / 64;
    delta = calloc(delta_words, sizeof(uint64_t));

    /* Process p=2 */
    process_prime_2(n);

    /* Process each relevant odd prime (p%8 in {5,7}) */
    for (uint64_t i = 1; i <= n / 2; i++) {
        if (is_prime[i]) {
            uint64_t p = 2 * i + 1;
            if ((p & 7) == 5 || (p & 7) == 7) {
                process_relevant_prime(p, n);
            }
        }
    }

    /* Final pass: compute running parity and sum factorials */
    uint64_t fact = 1;
    uint64_t ans = 0;
    int running_par = 0;
    for (uint64_t k = 1; k <= n; k++) {
        fact = (fact * k) % MOD;
        if (delta_get(k)) {
            running_par ^= 1;
        }
        if (running_par == 0) {
            ans += fact;
            if (ans >= MOD) ans -= MOD;
        }
    }

    free(is_prime);
    free(delta);
    return ans;
}

long long p937_native(void) {
    /* Self-tests from problem statement */
    if (compute_G(4) != 25) {
        fprintf(stderr, "self-test G(4) failed: got %llu\n",
                (unsigned long long)compute_G(4));
        return -1;
    }
    if (compute_G(7) != 745) {
        fprintf(stderr, "self-test G(7) failed\n");
        return -1;
    }
    uint64_t g100 = compute_G(100);
    if (g100 != 709772949ULL) {
        fprintf(stderr, "self-test G(100) failed: got %llu\n",
                (unsigned long long)g100);
        return -1;
    }

    return (long long)compute_G(100000000ULL);
}
