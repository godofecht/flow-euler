/* Project Euler 947: Fibonacci Residues
 * Compute S(10^6) modulo 999999893.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef uint32_t u32;
typedef int64_t i64;
typedef uint64_t u64;
typedef __int128 i128;

#define MOD 999999893LL
#define N 1000000

static u32 *spf;
static u32 *primes_arr;
static u32 num_primes_arr;

static i64 mod_pow(i64 a, i64 e, i64 mod) {
    i64 r = 1 % mod;
    a %= mod;
    if (a < 0) a += mod;
    while (e > 0) {
        if (e & 1) r = (i64)((__int128)r * a % mod);
        a = (i64)((__int128)a * a % mod);
        e >>= 1;
    }
    return r;
}

static void sieve_spf(u32 limit) {
    spf = (u32*)malloc((limit + 1) * sizeof(u32));
    for (u32 i = 0; i <= limit; i++) spf[i] = i;
    spf[0] = 0; spf[1] = 1;

    for (u32 i = 2; (u64)i * i <= limit; i++) {
        if (spf[i] == i) {
            for (u32 j = i * i; j <= limit; j += i) {
                if (spf[j] == j) spf[j] = i;
            }
        }
    }
}

/* Fast doubling: return (F_n, F_{n+1}) mod m */
static void fib_pair(u64 n, i64 mod, i64 *fn, i64 *fn1) {
    i64 a = 0, b = 1 % mod;
    if (n == 0) { *fn = 0; *fn1 = b; return; }

    int nbits = 64 - __builtin_clzll(n);
    for (int i = nbits - 1; i >= 0; i--) {
        i64 two_b_minus_a = (2 * b - a) % mod;
        if (two_b_minus_a < 0) two_b_minus_a += mod;
        i64 c = (i64)((__int128)a * two_b_minus_a % mod); /* F_{2k} */
        i64 d = (i64)((__int128)a * a % mod + (__int128)b * b % mod) % mod; /* F_{2k+1} */
        if (d < 0) d += mod;
        if ((n >> i) & 1) {
            a = d;
            b = (c + d) % mod;
        } else {
            a = c;
            b = d;
        }
    }
    *fn = a;
    *fn1 = b;
}

static int check_A_order(u64 n, i64 p) {
    i64 fn, fn1;
    fib_pair(n, p, &fn, &fn1);
    return fn == 0 && fn1 == 1;
}

static i64 pisano_prime(i64 p) {
    if (p == 2) return 3;
    if (p == 5) return 20;

    i64 residue = mod_pow(5, (p - 1) / 2, p);
    i64 candidate = (residue == 1) ? (p - 1) : (2 * (p + 1));

    /* Factor candidate using SPF */
    i64 x = candidate;
    i64 distinct_primes[32];
    int ndp = 0;
    while (x > 1) {
        i64 q = spf[x];
        distinct_primes[ndp++] = q;
        while (x % q == 0) x /= q;
    }

    i64 d = candidate;
    for (int i = 0; i < ndp; i++) {
        i64 q = distinct_primes[i];
        while (d % q == 0) {
            i64 nd = d / q;
            if (check_A_order((u64)nd, p)) {
                d = nd;
            } else {
                break;
            }
        }
    }
    return d;
}

static i64 has_short_period_factor(i64 p, i64 pi_p) {
    if (p == 5) return 5;
    if (p == 2) return 1;

    if (mod_pow(5, (p - 1) / 2, p) != 1) return 1;
    if (pi_p % 2 != 0) return 1;

    u64 n = (u64)(pi_p / 2);
    i64 fn, fn1;
    fib_pair(n, p, &fn, &fn1);
    i64 ln = (2 * fn1 - fn) % p;
    if (ln < 0) ln += p;
    i64 minus1_pow = (n % 2 == 0) ? 1 : (p - 1);
    i64 det = (1 + minus1_pow - ln) % p;
    if (det < 0) det += p;
    return (det == 0) ? 2 : 1;
}

/* Distribution entry: (period, count_mod_MOD) */
typedef struct {
    i64 period;
    i64 count;
} DistEntry;

/* Cache for prime power distributions */
typedef struct {
    i64 p, e;
    DistEntry entries[2];
    int num_entries;
} DistCacheEntry;

static DistCacheEntry *dist_cache = NULL;
static int dist_cache_size = 0;
static int dist_cache_cap = 0;

static DistEntry* get_distribution(i64 p, i64 e, i64 pi_p, i64 k, int *num_entries) {
    /* Check cache */
    for (int i = 0; i < dist_cache_size; i++) {
        if (dist_cache[i].p == p && dist_cache[i].e == e) {
            *num_entries = dist_cache[i].num_entries;
            return dist_cache[i].entries;
        }
    }

    /* Compute distribution */
    i64 pe1 = 1;
    for (i64 i = 0; i < e - 1; i++) pe1 *= p;
    i64 T = pi_p * pe1;
    i64 p2 = p * p;
    i64 total = 1;
    for (i64 i = 0; i < e - 1; i++) total *= p2;
    total *= (p2 - 1);

    DistEntry entries[2];
    int ne;

    if (k == 1) {
        entries[0].period = T;
        entries[0].count = total % MOD;
        ne = 1;
    } else if (k == 2) {
        i64 small_period = T / 2;
        i64 small_count = (p - 1) * pe1;
        i64 big_count = total - small_count;
        entries[0].period = small_period;
        entries[0].count = small_count % MOD;
        entries[1].period = T;
        entries[1].count = big_count % MOD;
        ne = 2;
    } else {
        /* k == 5 (only for p=5) */
        i64 small_period = T / 5;
        i64 small_count = (p - 1);
        for (i64 i = 0; i < e - 1; i++) small_count *= p2;
        i64 big_count = total - small_count;
        entries[0].period = small_period;
        entries[0].count = small_count % MOD;
        entries[1].period = T;
        entries[1].count = big_count % MOD;
        ne = 2;
    }

    /* Store in cache */
    if (dist_cache_size >= dist_cache_cap) {
        dist_cache_cap = dist_cache_cap ? dist_cache_cap * 2 : 1024;
        dist_cache = (DistCacheEntry*)realloc(dist_cache, dist_cache_cap * sizeof(DistCacheEntry));
    }
    dist_cache[dist_cache_size].p = p;
    dist_cache[dist_cache_size].e = e;
    dist_cache[dist_cache_size].num_entries = ne;
    for (int i = 0; i < ne; i++) dist_cache[dist_cache_size].entries[i] = entries[i];
    dist_cache_size++;

    *num_entries = ne;
    return dist_cache[dist_cache_size - 1].entries;
}

static i64 gcd_i64(i64 a, i64 b) {
    while (b) { i64 t = a % b; a = b; b = t; }
    return a;
}

/* Simple dict for combining distributions */
typedef struct {
    i64 key;
    i64 val;
} DictEntry;

static DictEntry cur_dict[256];
static DictEntry new_dict[256];

long long p947_native(void) {
    i64 max_spf = 2 * N + 2;
    sieve_spf((u32)max_spf);

    /* Precompute pi(p) and k for primes up to N */
    i64 *pi_prime = (i64*)calloc(N + 1, sizeof(i64));
    i64 *k_prime = (i64*)malloc((N + 1) * sizeof(i64));
    for (int i = 0; i <= N; i++) k_prime[i] = 1;

    for (i64 p = 2; p <= N; p++) {
        if (spf[p] == (u32)p) { /* p is prime */
            i64 pi_p = pisano_prime(p);
            pi_prime[p] = pi_p;
            k_prime[p] = has_short_period_factor(p, pi_p);
        }
    }

    i64 ans = 0;

    for (i64 n = 1; n <= N; n++) {
        i64 x = n;

        /* Collect distributions for each prime power */
        DistEntry *dists[8];
        int dist_counts[8];
        int num_dists = 0;

        while (x > 1) {
            i64 p = spf[x];
            i64 e = 0;
            while (x % p == 0) { x /= p; e++; }

            int ne;
            DistEntry *d = get_distribution(p, e, pi_prime[p], k_prime[p], &ne);
            dists[num_dists] = d;
            dist_counts[num_dists] = ne;
            num_dists++;
        }

        /* Combine distributions via CRT (lcm of periods, product of counts) */
        int cur_size = 1;
        cur_dict[0].key = 1;
        cur_dict[0].val = 1;

        for (int di = 0; di < num_dists; di++) {
            int new_size = 0;
            for (int ci = 0; ci < cur_size; ci++) {
                for (int ei = 0; ei < dist_counts[di]; ei++) {
                    i64 per1 = cur_dict[ci].key;
                    i64 per2 = dists[di][ei].period;
                    i64 g = gcd_i64(per1, per2);
                    i128 l = ((i128)per1 / g) * per2;
                    i64 v = (i64)((__int128)cur_dict[ci].val * dists[di][ei].count % MOD);

                    /* Find or insert l in new_dict */
                    int found = -1;
                    for (int j = 0; j < new_size; j++) {
                        if (new_dict[j].key == (i64)l) { found = j; break; }
                    }
                    if (found >= 0) {
                        new_dict[found].val = (new_dict[found].val + v) % MOD;
                    } else {
                        new_dict[new_size].key = (i64)l;
                        new_dict[new_size].val = v;
                        new_size++;
                    }
                }
            }
            /* Copy new_dict to cur_dict, removing zeros */
            cur_size = 0;
            for (int j = 0; j < new_size; j++) {
                if (new_dict[j].val != 0) {
                    cur_dict[cur_size++] = new_dict[j];
                }
            }
        }

        /* Compute Pn = sum of per^2 * count mod MOD */
        i64 Pn = 0;
        for (int j = 0; j < cur_size; j++) {
            i64 per2 = (i64)((__int128)cur_dict[j].key * cur_dict[j].key % MOD);
            Pn = (Pn + (i64)((__int128)per2 * cur_dict[j].val % MOD)) % MOD;
        }

        ans = (ans + (i64)((__int128)Pn * (N / n) % MOD)) % MOD;
    }

    free(spf);
    free(pi_prime);
    free(k_prime);
    free(dist_cache);

    return ans;
}
