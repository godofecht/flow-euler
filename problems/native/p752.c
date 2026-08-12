/* Project Euler 752: Powers of 1 + sqrt(7)
   Port of the Python reference solver. */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static int *spf = NULL;
static int *primes = NULL;
static int nprimes = 0;

static void smallest_prime_factors(int limit) {
    spf = (int*)calloc(limit + 1, sizeof(int));
    primes = (int*)malloc((limit / 2 + 100) * sizeof(int));
    nprimes = 0;

    for (int i = 2; i <= limit; i++) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes[nprimes++] = i;
        }
        for (int j = 0; j < nprimes; j++) {
            int p = primes[j];
            long long val = (long long)i * p;
            if (val > limit || p > spf[i]) break;
            spf[val] = p;
        }
    }
}

/* unique prime factors of n */
static int unique_prime_factors(int n, int *out) {
    int cnt = 0;
    while (n > 1) {
        int p = spf[n];
        out[cnt++] = p;
        while (n % p == 0) n /= p;
    }
    return cnt;
}

/* (1 + sqrt(7))^exponent mod modulus, returns (a, b) */
static void quadratic_power(long long exponent, long long modulus,
                            long long *ra, long long *rb) {
    long long res_a = 1, res_b = 0;
    long long base_a = 1, base_b = 1;

    while (exponent) {
        if (exponent & 1) {
            long long na = (base_a * res_a + 7 * base_b * res_b) % modulus;
            long long nb = (base_a * res_b + base_b * res_a) % modulus;
            res_a = na; res_b = nb;
        }
        long long na = (base_a * base_a + 7 * base_b * base_b) % modulus;
        long long nb = (2 * base_a * base_b) % modulus;
        base_a = na; base_b = nb;
        exponent >>= 1;
    }
    *ra = res_a; *rb = res_b;
}

static long long llgcd(long long a, long long b) {
    while (b) { long long t = a % b; a = b; b = t; }
    return a;
}

/* pow(7, (p-1)/2, p) */
static long long powmod(long long base, long long exp, long long mod) {
    long long r = 1 % mod;
    base %= mod;
    while (exp) {
        if (exp & 1) r = r * base % mod;
        base = base * base % mod;
        exp >>= 1;
    }
    return r;
}

static long long prime_order(int p) {
    if (p == 7) return 7;

    int factors[64];
    int nf;
    long long order;

    if (powmod(7, (p - 1) / 2, p) == 1) {
        order = p - 1;
        nf = unique_prime_factors(p - 1, factors);
    } else {
        order = (long long)(p - 1) * (p + 1);
        nf = unique_prime_factors(p - 1, factors);
        int extra[64];
        int ne = unique_prime_factors(p + 1, extra);
        for (int i = 0; i < ne; i++) factors[nf++] = extra[i];
    }

    /* deduplicate factors */
    for (int i = 0; i < nf; i++) {
        for (int j = i + 1; j < nf; j++) {
            if (factors[i] == factors[j]) { factors[j] = factors[nf-1]; nf--; j--; }
        }
    }

    for (int i = 0; i < nf; i++) {
        long long q = factors[i];
        while (order % q == 0) {
            long long ra, rb;
            quadratic_power(order / q, p, &ra, &rb);
            if (ra == 1 && rb == 0) {
                order /= q;
            } else {
                break;
            }
        }
    }
    return order;
}

static long long *orders = NULL;

static void build_prime_power_orders(int limit) {
    orders = (long long*)calloc(limit + 1, sizeof(long long));

    for (int pi = 0; pi < nprimes; pi++) {
        int p = primes[pi];
        if (p > limit) break;
        if (p < 5) continue;

        long long order = prime_order(p);
        long long pp = p;
        long long lifted = order;
        orders[pp] = lifted;

        while (pp * p <= limit) {
            pp *= p;
            long long ra, rb;
            quadratic_power(lifted, pp, &ra, &rb);
            if (!(ra == 1 && rb == 0)) {
                lifted *= p;
            }
            orders[pp] = lifted;
        }
    }
}

static long long G(int limit) {
    smallest_prime_factors(limit + 1);
    build_prime_power_orders(limit);

    long long *values = (long long*)calloc(limit + 1, sizeof(long long));
    long long total = 0;

    for (int n = 2; n <= limit; n++) {
        if (n % 2 == 0 || n % 3 == 0) continue;

        int p = spf[n];
        int remaining = n;
        long long pp = 1;
        while (remaining % p == 0) {
            remaining /= p;
            pp *= p;
        }

        long long prev = (remaining > 1) ? values[remaining] : 1;
        long long comp = orders[pp];
        long long g = llgcd(prev, comp);
        long long val = prev / g * comp;
        values[n] = val;
        total += val;
    }

    free(values);
    free(orders);
    free(spf);
    free(primes);
    return total;
}

long long p752_native(void) {
    return G(1000000);
}
