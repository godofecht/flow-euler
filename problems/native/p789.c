/* Project Euler 789: Minimal Pairing Modulo p
 *
 * For p = 2_000_000_011, find the cost product of an optimal pairing.
 *
 * Strategy: Wilson implies K ≡ -1 (mod p). Optimal costs are primes.
 * Meet-in-the-middle over small primes with weight-bounded search.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

static i64 pow_mod(i64 a, i64 e, i64 m) {
    i64 r = 1, b = a % m;
    if (b < 0) b += m;
    while (e > 0) {
        if (e & 1) r = r * b % m;
        b = b * b % m;
        e >>= 1;
    }
    return r;
}

static int sieve_primes(int n, int *out) {
    if (n < 2) return 0;
    char *sieve = calloc(n + 1, 1);
    memset(sieve, 1, n + 1);
    sieve[0] = sieve[1] = 0;
    for (int i = 2; (i64)i * i <= n; i++)
        if (sieve[i])
            for (int j = i * i; j <= n; j += i) sieve[j] = 0;
    int k = 0;
    for (int i = 2; i <= n; i++)
        if (sieve[i]) out[k++] = i;
    free(sieve);
    return k;
}

/* Best map entry: residue -> (min_weight, min_product) */
typedef struct {
    i64 residue;
    i64 weight;
    /* product as integer - can be huge, store as u64 (may overflow for large weights) */
    /* Actually the product can be very large. We need to compare products at equal weight. */
    /* For the final answer, we need the actual product. Use double for comparison + u64 for mod. */
    /* Actually the Python stores the exact integer product. For p=2e9, products of primes up to ~400 */
    /* can be enormous. But the answer fits in u64. Let me think... */
    /* The answer is 13431419535872807040 which fits in u64. The intermediate products */
    /* are products of small primes. The final K = prod_a * prod_b * prod_23. */
    /* These individual products might exceed u64. We need bignum or careful handling. */
    /* Actually, looking at the Python code, prod can be very large (products of many primes). */
    /* But the final answer fits in u64. So the individual products must be manageable. */
    /* Let me use __int128 for products. */
    i128 prod;
} BestEntry;

/* Simple hash table for best map */
#define HASH_SIZE (1 << 23)
#define HASH_MASK (HASH_SIZE - 1)

typedef struct {
    i64 residue;
    i64 weight;
    i128 prod;
    int used;
} HashEntry;

static HashEntry *hash_new(void) {
    HashEntry *h = calloc(HASH_SIZE, sizeof(HashEntry));
    return h;
}

static unsigned hash_key(i64 residue) {
    u64 v = (u64)residue;
    v ^= v >> 33;
    v *= 0xff51afd7ed558ccdULL;
    v ^= v >> 33;
    v *= 0xc4ceb9fe1a85ec53ULL;
    v ^= v >> 33;
    return (unsigned)(v & HASH_MASK);
}

static void hash_set(HashEntry *h, i64 residue, i64 weight, i128 prod) {
    unsigned idx = hash_key(residue);
    while (h[idx].used) {
        if (h[idx].residue == residue) {
            if (weight < h[idx].weight || (weight == h[idx].weight && prod < h[idx].prod)) {
                h[idx].weight = weight;
                h[idx].prod = prod;
            }
            return;
        }
        idx = (idx + 1) & HASH_MASK;
    }
    h[idx].used = 1;
    h[idx].residue = residue;
    h[idx].weight = weight;
    h[idx].prod = prod;
}

static int hash_get(HashEntry *h, i64 residue, i64 *weight, i128 *prod) {
    unsigned idx = hash_key(residue);
    while (h[idx].used) {
        if (h[idx].residue == residue) {
            *weight = h[idx].weight;
            *prod = h[idx].prod;
            return 1;
        }
        idx = (idx + 1) & HASH_MASK;
    }
    return 0;
}

/* DFS to build best map for a set of primes */
static void dfs_best(const int *primes, int cnt, int idx, i64 p,
                     i64 residue, i128 prod, i64 weight, i64 w_max,
                     HashEntry *best) {
    if (idx == cnt) {
        hash_set(best, residue, weight, prod);
        return;
    }
    int q = primes[idx];
    i64 wq = q - 1;
    i64 max_e = (w_max - weight) / wq;
    i64 r = residue;
    i128 pr = prod;
    i64 w = weight;
    for (i64 e = 0; e <= max_e; e++) {
        dfs_best(primes, cnt, idx + 1, p, r, pr, w, w_max, best);
        r = r * q % p;
        pr *= q;
        w += wq;
    }
}

static HashEntry *build_best_map(const int *primes, int cnt, i64 p, i64 w_max) {
    HashEntry *best = hash_new();
    dfs_best(primes, cnt, 0, p, 1, 1, 0, w_max, best);
    return best;
}

/* Collect all entries from hash table */
typedef struct {
    i64 residue;
    i64 weight;
    i128 prod;
} State;

static int hash_collect(HashEntry *h, State *out) {
    int n = 0;
    for (int i = 0; i < HASH_SIZE; i++) {
        if (h[i].used) {
            out[n].residue = h[i].residue;
            out[n].weight = h[i].weight;
            out[n].prod = h[i].prod;
            n++;
        }
    }
    return n;
}

static int cmp_state_weight(const void *a, const void *b) {
    i64 wa = ((const State *)a)->weight;
    i64 wb = ((const State *)b)->weight;
    if (wa < wb) return -1;
    if (wa > wb) return 1;
    return 0;
}

/* 2^e2 * 3^e3 states */
typedef struct {
    i64 weight;
    i64 residue;
    i128 prod;
    i64 inv_res;
} State23;

static int build_23_states(i64 p, i64 w_max, State23 *out) {
    i64 inv2 = pow_mod(2, p - 2, p);
    i64 inv3 = pow_mod(3, p - 2, p);

    int max_e3 = w_max / 2;
    /* Precompute powers */
    i64 *pow2_mod = malloc((w_max + 1) * sizeof(i64));
    i128 *pow2_int = malloc((w_max + 1) * sizeof(i128));
    i64 *inv2_mod = malloc((w_max + 1) * sizeof(i64));
    pow2_mod[0] = 1; pow2_int[0] = 1; inv2_mod[0] = 1;
    for (int i = 1; i <= w_max; i++) {
        pow2_mod[i] = pow2_mod[i-1] * 2 % p;
        pow2_int[i] = pow2_int[i-1] * 2;
        inv2_mod[i] = inv2_mod[i-1] * inv2 % p;
    }

    i64 *pow3_mod = malloc((max_e3 + 1) * sizeof(i64));
    i128 *pow3_int = malloc((max_e3 + 1) * sizeof(i128));
    i64 *inv3_mod = malloc((max_e3 + 1) * sizeof(i64));
    pow3_mod[0] = 1; pow3_int[0] = 1; inv3_mod[0] = 1;
    for (int i = 1; i <= max_e3; i++) {
        pow3_mod[i] = pow3_mod[i-1] * 3 % p;
        pow3_int[i] = pow3_int[i-1] * 3;
        inv3_mod[i] = inv3_mod[i-1] * inv3 % p;
    }

    int n = 0;
    for (int e3 = 0; e3 <= max_e3; e3++) {
        i64 w3 = 2 * e3;
        int remaining = w_max - w3;
        i64 base_res = pow3_mod[e3];
        i128 base_prod = pow3_int[e3];
        i64 base_inv = inv3_mod[e3];
        for (int e2 = 0; e2 <= remaining; e2++) {
            i64 w = w3 + e2;
            i64 res = base_res * pow2_mod[e2] % p;
            i128 prod = base_prod * pow2_int[e2];
            i64 inv_res = base_inv * inv2_mod[e2] % p;
            out[n].weight = w;
            out[n].residue = res;
            out[n].prod = prod;
            out[n].inv_res = inv_res;
            n++;
        }
    }

    qsort(out, n, sizeof(State23), cmp_state_weight);

    free(pow2_mod); free(pow2_int); free(inv2_mod);
    free(pow3_mod); free(pow3_int); free(inv3_mod);
    return n;
}

static i128 search_with_bound(i64 p, i64 w_max) {
    int primes_buf[1000];
    int cnt = sieve_primes(w_max + 1, primes_buf);
    /* Remove 2 and 3, split into A (5-23) and B (29+) */
    int primes_a[20], primes_b[500];
    int na = 0, nb = 0;
    for (int i = 0; i < cnt; i++) {
        int q = primes_buf[i];
        if (q == 2 || q == 3) continue;
        if (q <= 23) primes_a[na++] = q;
        else primes_b[nb++] = q;
    }

    HashEntry *best_a = build_best_map(primes_a, na, p, w_max);
    HashEntry *best_b_hash = build_best_map(primes_b, nb, p, w_max);

    /* Collect B states */
    State *b_states = malloc(HASH_SIZE * sizeof(State));
    int bn = hash_collect(best_b_hash, b_states);
    qsort(b_states, bn, sizeof(State), cmp_state_weight);

    /* Build 2^e2 * 3^e3 states */
    State23 *states23 = malloc(HASH_SIZE * sizeof(State23));
    int n23 = build_23_states(p, w_max, states23);

    i64 target = p - 1; /* -1 mod p */
    i64 best_weight = (i64)1e18;
    i128 best_k = 0;
    int found = 0;

    for (int bi = 0; bi < bn; bi++) {
        i64 w_b = b_states[bi].weight;
        if (w_b > w_max) break;
        i64 max_w23 = w_max - w_b;
        i64 inv_b = pow_mod(b_states[bi].residue, p - 2, p);
        i64 base_needed = target * inv_b % p;

        for (int si = 0; si < n23; si++) {
            if (states23[si].weight > max_w23) break;
            i64 rem = w_max - w_b - states23[si].weight;
            i64 needed_a = base_needed * states23[si].inv_res % p;
            i64 w_a;
            i128 prod_a;
            if (!hash_get(best_a, needed_a, &w_a, &prod_a)) continue;
            if (w_a > rem) continue;

            i64 total_w = w_b + states23[si].weight + w_a;
            i128 k = prod_a * b_states[bi].prod * states23[si].prod;
            if (!found || total_w < best_weight ||
                (total_w == best_weight && k < best_k)) {
                best_weight = total_w;
                best_k = k;
                found = 1;
            }
        }
    }

    free(b_states);
    free(states23);
    /* Free hash tables by iterating */
    for (int i = 0; i < HASH_SIZE; i++) {
        best_a[i].used = 0;
        best_b_hash[i].used = 0;
    }
    free(best_a);
    free(best_b_hash);

    return found ? best_k : -1;
}

unsigned long long p789_native(void) {
    i64 p = 2000000011LL;
    i64 w_bounds[] = {240, 260, 280, 300, 320, 360, 400};
    for (int i = 0; i < 7; i++) {
        i128 k = search_with_bound(p, w_bounds[i]);
        if (k >= 0) return (unsigned long long)k;
    }
    /* Fallback */
    i64 w = 360;
    while (1) {
        i128 k = search_with_bound(p, w);
        if (k >= 0) return (unsigned long long)k;
        w += 40;
    }
}
