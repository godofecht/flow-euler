/* Project Euler 753: Fermat Equation
   Sum over primes p < 6,000,000 of F(p),
   where F(p) counts ordered triples (a,b,c) with 1<=a,b,c<p and a^3+b^3≡c^3 (mod p).
   Port of the Python reference solver. */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

static int isqrt_ll(long long n) {
    if (n <= 0) return 0;
    long long r = (long long)sqrt((double)n);
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return (int)r;
}

/* Odd-only sieve up to n (inclusive).
   sieve_odd[i] is 1 iff (2*i+1) is prime (for i>=1).
   sieve_odd[0] corresponds to 1 (not prime). */
static unsigned char *sieve_odd = NULL;
static int *primes_list = NULL;
static int nprimes = 0;

static void sieve_primes_upto(int n) {
    if (n < 2) { sieve_odd = NULL; primes_list = NULL; nprimes = 0; return; }

    int size = (n / 2) + 1;
    sieve_odd = (unsigned char*)malloc(size);
    memset(sieve_odd, 1, size);
    sieve_odd[0] = 0; /* 1 is not prime */

    int r = isqrt_ll(n);
    for (int i = 1; i <= r / 2; i++) {
        if (sieve_odd[i]) {
            int p = 2 * i + 1;
            int start = (p * p) / 2;
            for (int j = start; j < size; j += p)
                sieve_odd[j] = 0;
        }
    }

    /* count primes */
    int cnt = 1; /* 2 is prime */
    for (int i = 1; i < size; i++) if (sieve_odd[i]) cnt++;

    primes_list = (int*)malloc(cnt * sizeof(int));
    nprimes = 0;
    primes_list[nprimes++] = 2;
    for (int i = 1; i < size; i++) {
        if (sieve_odd[i]) primes_list[nprimes++] = 2 * i + 1;
    }
    /* remove last if > n */
    if (nprimes > 0 && primes_list[nprimes - 1] > n) nprimes--;
}

/* Build u_by_p for primes p ≡ 1 (mod 3) with 4p = u^2 + 27 v^2 */
static unsigned short *u_by_p = NULL;

static void build_u_map(int limit) {
    int n = limit;
    long long n4 = 4LL * n;
    int vmax = isqrt_ll(n4 / 27);
    u_by_p = (unsigned short*)calloc(n + 1, sizeof(unsigned short));

    for (int v = 1; v <= vmax; v++) {
        long long base = 27LL * v * v;
        int umax = isqrt_ll(n4 - base);
        /* u parity must match v parity */
        int u = v & 1;
        while (u <= umax) {
            long long p = (u * u + base) >> 2;
            if (p <= n && (p % 3 == 1)) {
                if (p != 2 && (p & 1) && sieve_odd[p >> 1]) {
                    if (u_by_p[p] == 0) u_by_p[p] = (unsigned short)u;
                }
            }
            u += 2;
        }
    }
}

static int trace_ap_for_curve(int p) {
    /* p ≡ 1 (mod 3), p != 3 */
    int u = u_by_p[p];
    /* u should be nonzero */
    return (u % 3 == 2) ? u : -u;
}

static long long F_of_prime(int p) {
    if (p == 3) {
        /* direct count for p=3 */
        int cubes[3] = {0, 1, 2};
        int cnt = 0;
        for (int a = 1; a <= 2; a++) {
            for (int b = 1; b <= 2; b++) {
                int s = (cubes[a] + cubes[b]) % 3;
                for (int c = 1; c <= 2; c++) {
                    if (cubes[c] == s) cnt++;
                }
            }
        }
        return cnt;
    }

    if (p % 3 == 2) {
        return (long long)(p - 1) * (p - 2);
    }

    /* p ≡ 1 (mod 3) */
    int ap = trace_ap_for_curve(p);
    return (long long)(p - 1) * (p - ap - 8);
}

long long p753_native(void) {
    int limit_exclusive = 6000000;
    int max_p = limit_exclusive - 1;

    sieve_primes_upto(max_p);
    build_u_map(max_p);

    long long total = 0;
    for (int i = 0; i < nprimes; i++) {
        int p = primes_list[i];
        if (p >= limit_exclusive) break;
        total += F_of_prime(p);
    }

    free(sieve_odd);
    free(primes_list);
    free(u_by_p);
    return total;
}
