// Project Euler 758: Buckets of Water.
//
// Sum over prime pairs p<q<1000 of P(2^{p^5}-1, 2^{q^5}-1) mod 1e9+7.
// P(a,b) is the minimal pour count to measure 1 litre with buckets
// (a, b, a+b). A continued-fraction characterization gives
//   P(a,b) = 2*(pen_p + pen_q) - 2
// where pen_p/pen_q is the penultimate convergent of b/a.
//
// For Mersenne arguments (2^E-1)/(2^F-1) the Euclidean step reduces to
// Euclid on the exponents, so the CF terms are computed mod MOD without
// ever forming the huge integers.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static const long long MOD = 1000000007LL;

static long long mulmod(long long a, long long b) {
    return (a % MOD) * (b % MOD) % MOD;
}

static long long powmod(long long base, long long exp) {
    long long result = 1 % MOD;
    base %= MOD;
    if (base < 0) base += MOD;
    while (exp > 0) {
        if (exp & 1) result = mulmod(result, base);
        base = mulmod(base, base);
        exp >>= 1;
    }
    return result;
}

static long long invmod(long long a) {
    // Fermat: a^(MOD-2) mod MOD
    return powmod(a % MOD, MOD - 2);
}

// S = 1 + ratio + ... + ratio^(m-1) (mod MOD), m >= 1.
static long long geom_sum_ratio(long long ratio, long long m) {
    if (m <= 0) return 0;
    ratio %= MOD;
    if (ratio < 0) ratio += MOD;
    if (ratio == 1) return m % MOD;
    long long rm = powmod(ratio, m);
    long long num = (rm - 1) % MOD;
    if (num < 0) num += MOD;
    return mulmod(num, invmod(ratio - 1));
}

// Continued-fraction terms of (2^e_large - 1)/(2^e_small - 1), each mod MOD.
// Writes into terms[] and returns the count.
static int cf_terms_mersenne(long long e_small, long long e_large,
                             long long *terms, int cap) {
    int n = 0;
    long long hi = e_large, lo = e_small;
    for (;;) {
        long long m = hi / lo;
        long long r = hi % lo;
        long long ratio = powmod(2, lo);
        long long shift = powmod(2, r);
        long long series = geom_sum_ratio(ratio, m);
        long long q_mod = mulmod(shift, series);
        terms[n++] = q_mod;
        if (r == 0) break;
        hi = lo;
        lo = r;
    }
    return n;
}

// Penultimate convergent (p, q) of the CF with the given terms, mod MOD.
// If only one term, returns that convergent.
static void penultimate_convergent(const long long *terms, int n,
                                   long long *out_p, long long *out_q) {
    long long p_m2 = 0, p_m1 = 1;
    long long q_m2 = 1, q_m1 = 0;
    long long last_p = 0, last_q = 0, prev_p = 0, prev_q = 0;
    for (int i = 0; i < n; i++) {
        long long a = terms[i];
        long long p = (mulmod(a, p_m1) + p_m2) % MOD;
        long long q = (mulmod(a, q_m1) + q_m2) % MOD;
        prev_p = last_p; prev_q = last_q;
        last_p = p; last_q = q;
        p_m2 = p_m1; p_m1 = p;
        q_m2 = q_m1; q_m1 = q;
    }
    if (n == 1) {
        *out_p = last_p; *out_q = last_q;
    } else {
        *out_p = prev_p; *out_q = prev_q;
    }
}

static long long P_mersenne_exponents(long long ea, long long eb) {
    long long terms[256];
    int n = cf_terms_mersenne(ea, eb, terms, 256);
    long long p, q;
    penultimate_convergent(terms, n, &p, &q);
    long long s = (p + q) % MOD;
    long long result = (2 * s - 2) % MOD;
    if (result < 0) result += MOD;
    return result;
}

static int primes_below(int n, int *out) {
    // sieve
    char *sieve = (char *)calloc(n, 1);
    for (int i = 2; i < n; i++) sieve[i] = 1;
    for (int p = 2; (long long)p * p < n; p++) {
        if (sieve[p]) {
            for (int k = p * p; k < n; k += p) sieve[k] = 0;
        }
    }
    int count = 0;
    for (int i = 2; i < n; i++) if (sieve[i]) out[count++] = i;
    free(sieve);
    return count;
}

long long p758_native(void) {
    int primes[1000];
    int np = primes_below(1000, primes);

    long long exps[1000];
    for (int i = 0; i < np; i++) {
        long long p = primes[i];
        long long e = p;
        // p^5
        e = e * p * p * p * p;  // p * p^4 = p^5; fits since p<1000 -> p^5 < 1e15
        exps[i] = e;
    }

    long long total = 0;
    for (int i = 0; i < np; i++) {
        for (int j = i + 1; j < np; j++) {
            total = (total + P_mersenne_exponents(exps[i], exps[j])) % MOD;
        }
    }
    return total;
}
