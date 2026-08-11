// Project Euler 787: Bezout's Game
// Port of the Python reference solver to C.
// Counts winning positions (a,b) with gcd(a,b)=1, a>0, b>0, a+b<=N.
//
// Derived rule:
//   a+b even  -> always winning
//   a+b odd   -> losing iff min(a,b) is even
//
// H(N) = total_coprime_pairs - losing_pairs
// Total coprime ordered pairs with sum<=N is sum_{s=2..N} phi(s) = S_phi(N)-1.
// Losing unordered pairs counted via Mobius inversion over odd d.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef long long i64;
typedef unsigned long long u64;

/* ---- Open-addressing hash map for Du Jiao sieve memoization ---- */

#define HASH_SIZE (1 << 18)
#define HASH_MASK (HASH_SIZE - 1)

typedef struct {
    u64 key;
    i64 val;
} HashEntry;

typedef struct {
    HashEntry *entries;
    char *used;
} HashMap;

static HashMap *hm_new(void) {
    HashMap *hm = (HashMap *)malloc(sizeof(HashMap));
    hm->entries = (HashEntry *)calloc(HASH_SIZE, sizeof(HashEntry));
    hm->used = (char *)calloc(HASH_SIZE, 1);
    return hm;
}

static u64 hash64(u64 x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

static i64 hm_get(HashMap *hm, i64 key, int *found) {
    u64 h = hash64((u64)key) & HASH_MASK;
    while (hm->used[h]) {
        if (hm->entries[h].key == (u64)key) {
            *found = 1;
            return hm->entries[h].val;
        }
        h = (h + 1) & HASH_MASK;
    }
    *found = 0;
    return 0;
}

static void hm_put(HashMap *hm, i64 key, i64 val) {
    u64 h = hash64((u64)key) & HASH_MASK;
    while (hm->used[h]) {
        if (hm->entries[h].key == (u64)key) {
            hm->entries[h].val = val;
            return;
        }
        h = (h + 1) & HASH_MASK;
    }
    hm->entries[h].key = (u64)key;
    hm->entries[h].val = val;
    hm->used[h] = 1;
}

/* ---- Linear sieve: mu, phi, and prefix sums ---- */

static i64 L;
static signed char *mu;
static int *phi;
static i64 *pre_mu;
static i64 *pre_phi;
static i64 *pre_mu_odd;

static HashMap *memo_mu;
static HashMap *memo_phi;
static HashMap *memo_mu_odd;

static void sieve(i64 limit) {
    mu = (signed char *)calloc((size_t)(limit + 1), 1);
    phi = (int *)calloc((size_t)(limit + 1), sizeof(int));
    char *is_comp = (char *)calloc((size_t)(limit + 1), 1);
    i64 *primes = (i64 *)malloc(sizeof(i64) * (size_t)limit);
    i64 np = 0;

    mu[1] = 1;
    phi[1] = 1;

    for (i64 i = 2; i <= limit; i++) {
        if (!is_comp[i]) {
            primes[np++] = i;
            mu[i] = -1;
            phi[i] = (int)(i - 1);
        }
        for (i64 j = 0; j < np; j++) {
            i64 p = primes[j];
            i64 v = i * p;
            if (v > limit) break;
            is_comp[v] = 1;
            if (i % p == 0) {
                mu[v] = 0;
                phi[v] = phi[i] * (int)p;
                break;
            } else {
                mu[v] = -mu[i];
                phi[v] = phi[i] * (int)(p - 1);
            }
        }
    }

    pre_mu = (i64 *)calloc((size_t)(limit + 1), sizeof(i64));
    pre_phi = (i64 *)calloc((size_t)(limit + 1), sizeof(i64));
    pre_mu_odd = (i64 *)calloc((size_t)(limit + 1), sizeof(i64));

    for (i64 i = 1; i <= limit; i++) {
        pre_mu[i] = pre_mu[i - 1] + mu[i];
        pre_phi[i] = pre_phi[i - 1] + phi[i];
        pre_mu_odd[i] = pre_mu_odd[i - 1] + ((i & 1) ? mu[i] : 0);
    }

    free(is_comp);
    free(primes);
}

/* ---- Du Jiao sieve for summatory functions ---- */

static i64 S_mu(i64 n) {
    if (n <= L) return pre_mu[n];
    int found;
    i64 v = hm_get(memo_mu, n, &found);
    if (found) return v;
    i64 res = 1;
    i64 i = 2;
    while (i <= n) {
        i64 q = n / i;
        i64 j = n / q;
        res -= (j - i + 1) * S_mu(q);
        i = j + 1;
    }
    hm_put(memo_mu, n, res);
    return res;
}

static i64 S_phi(i64 n) {
    if (n <= L) return pre_phi[n];
    int found;
    i64 v = hm_get(memo_phi, n, &found);
    if (found) return v;
    i64 res = n * (n + 1) / 2;
    i64 i = 2;
    while (i <= n) {
        i64 q = n / i;
        i64 j = n / q;
        res -= (j - i + 1) * S_phi(q);
        i = j + 1;
    }
    hm_put(memo_phi, n, res);
    return res;
}

static i64 S_mu_odd(i64 n) {
    if (n <= L) return pre_mu_odd[n];
    int found;
    i64 v = hm_get(memo_mu_odd, n, &found);
    if (found) return v;
    /* M_odd(n) = M(n) + M_odd(floor(n/2)) */
    i64 res = S_mu(n) + S_mu_odd(n / 2);
    hm_put(memo_mu_odd, n, res);
    return res;
}

/* ---- Counting helpers ---- */

/* Count (x,y) with x>=1, y odd>=1, 2x<y, 2x+y<=M. */
static i64 C_count(i64 M) {
    i64 t = (M - 1) / 4;
    if (t <= 0) return 0;
    return t * ((M + 1) / 2) - t * (t + 1);
}

static i64 H(i64 N) {
    L = (i64)(pow((double)N, 2.0 / 3.0)) + 10;
    if (L > N) L = N;

    sieve(L);

    memo_mu = hm_new();
    memo_phi = hm_new();
    memo_mu_odd = hm_new();

    /* Force-fill S_mu memo for N (speeds up S_mu_odd queries). */
    S_mu(N);

    /* Total coprime ordered pairs with sum<=N is sum_{s=2..N} phi(s). */
    i64 total_positions = S_phi(N) - 1;

    /* Losing unordered pairs: sum_{d odd} mu(d) * C(floor(N/d)). */
    i64 losing_unordered = 0;
    i64 d = 1;
    while (d <= N) {
        i64 q = N / d;
        i64 nd = N / q;
        i64 mu_range_odd = S_mu_odd(nd) - S_mu_odd(d - 1);
        losing_unordered += mu_range_odd * C_count(q);
        d = nd + 1;
    }

    i64 losing_ordered = 2 * losing_unordered;
    return total_positions - losing_ordered;
}

long long p787_native(void) {
    return H(1000000000LL);
}
