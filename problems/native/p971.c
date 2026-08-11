/* Project Euler 971
   For each prime p with p % 5 == 1, let n = (p-1)/5.  Consider the map
   phi(s) = s * (1+s)^n mod p on the 5-element subgroup mu5 of F_p^*.
   t = number of elements of mu5 lying on a cycle of phi.
   C_p = 1 + n * t.  S = sum of C_p over primes p <= 10^8 with p % 5 == 1.
   Answer: 33626723890930.
*/

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

typedef int64_t i64;
typedef uint64_t u64;
typedef __int128 i128;

static u64 mulmod(u64 a, u64 b, u64 m) {
    return (u64)((__uint128_t)a * b % m);
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

/* Odd-only sieve up to n. Returns array of primes and count. */
static int sieve_primes(u64 n, u64 **out_primes) {
    if (n < 2) { *out_primes = NULL; return 0; }
    u64 half = n / 2 + 1;          /* indices 0..half-1 represent 2*i+1 */
    unsigned char *sv = (unsigned char *)calloc(half, 1);
    if (!sv) { *out_primes = NULL; return 0; }
    u64 limit = (u64)sqrtl((long double)n);
    for (u64 i = 1; i <= limit / 2; i++) {
        if (sv[i] == 0) {
            u64 p = 2 * i + 1;
            u64 start = (p * p - 1) / 2;
            for (u64 j = start; j < half; j += p) sv[j] = 1;
        }
    }
    /* count */
    u64 cnt = 1; /* prime 2 */
    for (u64 i = 1; i < half; i++) {
        u64 v = 2 * i + 1;
        if (v <= n && sv[i] == 0) cnt++;
    }
    u64 *primes = (u64 *)malloc(cnt * sizeof(u64));
    if (!primes) { free(sv); *out_primes = NULL; return 0; }
    u64 k = 0;
    primes[k++] = 2;
    for (u64 i = 1; i < half; i++) {
        u64 v = 2 * i + 1;
        if (v <= n && sv[i] == 0) primes[k++] = v;
    }
    free(sv);
    *out_primes = primes;
    return (int)k;
}

/* For prime p with p % 5 == 1, count elements of mu5 on a cycle of phi. */
static int count_periodic_in_mu5(u64 p) {
    u64 n = (p - 1) / 5;
    /* find a generator of mu5: a with a^n != 1 mod p */
    u64 a = 2;
    u64 zeta = powmod(a, n, p);
    while (zeta == 1) {
        a++;
        zeta = powmod(a, n, p);
    }
    /* build mu5: {1, zeta, zeta^2, zeta^3, zeta^4} */
    u64 mu5[5];
    mu5[0] = 1 % p;
    u64 cur = 1;
    for (int i = 1; i < 5; i++) {
        cur = mulmod(cur, zeta, p);
        mu5[i] = cur;
    }
    /* phi(s) = s * (1+s)^n mod p */
    u64 phi[5];
    for (int i = 0; i < 5; i++) {
        u64 s = mu5[i];
        u64 one_plus = (1 + s) % p;
        u64 pv = powmod(one_plus, n, p);
        phi[i] = mulmod(s, pv, p);
    }
    /* map each mu5 element to its index; others map to -1 */
    /* functional graph traversal on the 5 nodes */
    int visited[5];        /* 0 unvisited, 1 in stack, 2 done */
    int in_cycle[5];
    for (int i = 0; i < 5; i++) { visited[i] = 0; in_cycle[i] = 0; }

    /* helper: find index of value v in mu5, or -1 */
    for (int start = 0; start < 5; start++) {
        if (visited[start]) continue;
        int path[8];
        int pathlen = 0;
        int ci = start;
        while (1) {
            if (visited[ci] == 0) {
                visited[ci] = 1;
                path[pathlen++] = ci;
                /* compute next: phi[ci]; find its index */
                u64 v = phi[ci];
                int next = -1;
                for (int j = 0; j < 5; j++) {
                    if (mu5[j] == v) { next = j; break; }
                }
                if (next < 0) {
                    /* maps outside mu5: no cycle on this path */
                    for (int x = 0; x < pathlen; x++) visited[path[x]] = 2;
                    break;
                }
                ci = next;
            } else if (visited[ci] == 1) {
                /* found a cycle: from first occurrence of ci in path */
                int idx = 0;
                for (int x = 0; x < pathlen; x++) {
                    if (path[x] == ci) { idx = x; break; }
                }
                for (int x = idx; x < pathlen; x++) in_cycle[path[x]] = 1;
                for (int x = 0; x < pathlen; x++) visited[path[x]] = 2;
                break;
            } else {
                /* visited == 2: leads to processed region */
                for (int x = 0; x < pathlen; x++) visited[path[x]] = 2;
                break;
            }
        }
    }
    int t = 0;
    for (int i = 0; i < 5; i++) t += in_cycle[i];
    return t;
}

long long p971_native(void) {
    u64 limit = 100000000ULL; /* 10^8 */
    u64 *primes;
    int np = sieve_primes(limit, &primes);
    i128 total = 0;
    for (int i = 0; i < np; i++) {
        u64 p = primes[i];
        if (p < 5) {
            if (p == 2 || p == 3) continue;
        }
        if (p % 5 != 1) continue;
        u64 n = (p - 1) / 5;
        int t = count_periodic_in_mu5(p);
        i64 C_p = 1 + (i64)n * t;
        total += C_p;
    }
    free(primes);
    return (long long)total;
}
