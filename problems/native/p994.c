// Project Euler 994: Counting Triangles
// T(1234*10^8, 2345*10^8) mod 1e9+7.
// Port of the Python reference solver to C.
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef int64_t i64;
typedef __int128 i128;

#define MOD 1000000007LL

static i64 INV2, INV6;

static i64 mulmod(i64 a, i64 b) {
    i128 r = (i128)a * b % MOD;
    return (i64)(r < 0 ? r + MOD : r);
}

static i64 mod_pow(i64 a, i64 e) {
    i64 r = 1;
    a %= MOD;
    if (a < 0) a += MOD;
    while (e > 0) {
        if (e & 1) r = mulmod(r, a);
        a = mulmod(a, a);
        e >>= 1;
    }
    return r;
}

static i64 p1_func(i64 n) {
    i64 a = n % MOD;
    i64 b = (a + 1) % MOD;
    return mulmod(mulmod(a, b), INV2);
}

static i64 p2_func(i64 n) {
    i64 a = n % MOD;
    i64 b = (a + 1) % MOD;
    i64 c = (2 * a + 1) % MOD;
    return mulmod(mulmod(mulmod(a, b), c), INV6);
}

static i64 p3_func(i64 n) {
    i64 s = p1_func(n);
    return mulmod(s, s);
}

static i64 c2_mod(i64 x) {
    i64 a = x % MOD;
    i64 b = (a - 1 + MOD) % MOD;
    return mulmod(mulmod(a, b), INV2);
}

static i64 c3_mod(i64 x) {
    i64 a = x % MOD;
    i64 b = (a - 1 + MOD) % MOD;
    i64 c = (a - 2 + MOD) % MOD;
    return mulmod(mulmod(mulmod(a, b), c), INV6);
}

/* ---- Hash table for Du Jiao cache ---- */
typedef struct {
    i64 key;
    i64 f0, f1, f2;
} CacheEntry;

typedef struct {
    CacheEntry *entries;
    uint8_t *used;
    i64 size;
    i64 count;
} Cache;

static i64 cache_hash(i64 key, i64 size) {
    uint64_t k = (uint64_t)key;
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;
    return (i64)(k % (uint64_t)size);
}

static void cache_put(Cache *c, i64 key, i64 f0, i64 f1, i64 f2);

static void cache_resize(Cache *c, i64 new_size) {
    CacheEntry *old_e = c->entries;
    uint8_t *old_u = c->used;
    i64 old_s = c->size;
    c->entries = (CacheEntry *)calloc((size_t)new_size, sizeof(CacheEntry));
    c->used = (uint8_t *)calloc((size_t)new_size, 1);
    c->size = new_size;
    c->count = 0;
    for (i64 i = 0; i < old_s; i++)
        if (old_u[i])
            cache_put(c, old_e[i].key, old_e[i].f0, old_e[i].f1, old_e[i].f2);
    free(old_e);
    free(old_u);
}

static void cache_put(Cache *c, i64 key, i64 f0, i64 f1, i64 f2) {
    if (c->count * 2 >= c->size)
        cache_resize(c, c->size * 2);
    i64 idx = cache_hash(key, c->size);
    while (c->used[idx]) {
        if (c->entries[idx].key == key) {
            c->entries[idx].f0 = f0;
            c->entries[idx].f1 = f1;
            c->entries[idx].f2 = f2;
            return;
        }
        idx = (idx + 1) % c->size;
    }
    c->entries[idx].key = key;
    c->entries[idx].f0 = f0;
    c->entries[idx].f1 = f1;
    c->entries[idx].f2 = f2;
    c->used[idx] = 1;
    c->count++;
}

static int cache_get(Cache *c, i64 key, i64 *f0, i64 *f1, i64 *f2) {
    i64 idx = cache_hash(key, c->size);
    while (c->used[idx]) {
        if (c->entries[idx].key == key) {
            *f0 = c->entries[idx].f0;
            *f1 = c->entries[idx].f1;
            *f2 = c->entries[idx].f2;
            return 1;
        }
        idx = (idx + 1) % c->size;
    }
    return 0;
}

/* ---- Totient prefix sums with Du Jiao recursion ---- */
typedef struct {
    i64 limit;
    uint32_t *pref0;
    uint32_t *pref1;
    uint32_t *pref2;
    Cache cache;
} TotientPrefix;

static void tp_values(TotientPrefix *tp, i64 n,
                       i64 *f0, i64 *f1, i64 *f2) {
    if (n <= tp->limit) {
        *f0 = tp->pref0[n];
        *f1 = tp->pref1[n];
        *f2 = tp->pref2[n];
        return;
    }
    if (cache_get(&tp->cache, n, f0, f1, f2)) return;

    *f0 = p1_func(n);
    *f1 = p2_func(n);
    *f2 = p3_func(n);

    i64 l = 2;
    while (l <= n) {
        i64 q = n / l;
        i64 r = n / q;

        i64 sum0 = (r - l + 1) % MOD;
        i64 sum1 = (p1_func(r) - p1_func(l - 1) + MOD) % MOD;
        i64 sum2 = (p2_func(r) - p2_func(l - 1) + MOD) % MOD;

        i64 sub0, sub1, sub2;
        tp_values(tp, q, &sub0, &sub1, &sub2);

        *f0 = (*f0 - mulmod(sum0, sub0) + MOD) % MOD;
        *f1 = (*f1 - mulmod(sum1, sub1) + MOD) % MOD;
        *f2 = (*f2 - mulmod(sum2, sub2) + MOD) % MOD;

        l = r + 1;
    }

    cache_put(&tp->cache, n, *f0, *f1, *f2);
}

static void tp_init(TotientPrefix *tp, i64 limit) {
    tp->limit = limit;

    uint32_t *phi = (uint32_t *)malloc((size_t)(limit + 1) * sizeof(uint32_t));
    for (i64 i = 0; i <= limit; i++) phi[i] = (uint32_t)i;
    for (i64 p = 2; p <= limit; p++) {
        if (phi[p] == (uint32_t)p) {
            for (i64 j = p; j <= limit; j += p)
                phi[j] -= phi[j] / (uint32_t)p;
        }
    }

    tp->pref0 = (uint32_t *)calloc((size_t)(limit + 1), sizeof(uint32_t));
    tp->pref1 = (uint32_t *)calloc((size_t)(limit + 1), sizeof(uint32_t));
    tp->pref2 = (uint32_t *)calloc((size_t)(limit + 1), sizeof(uint32_t));

    i64 s0 = 0, s1 = 0, s2 = 0;
    for (i64 i = 1; i <= limit; i++) {
        i64 ph = phi[i];
        i64 im = i % MOD;
        s0 = (s0 + ph) % MOD;
        s1 = (s1 + im * ph % MOD) % MOD;
        s2 = (s2 + im * im % MOD * ph % MOD) % MOD;
        tp->pref0[i] = (uint32_t)s0;
        tp->pref1[i] = (uint32_t)s1;
        tp->pref2[i] = (uint32_t)s2;
    }
    free(phi);

    tp->cache.size = 1 << 20;
    tp->cache.count = 0;
    tp->cache.entries = (CacheEntry *)calloc((size_t)tp->cache.size, sizeof(CacheEntry));
    tp->cache.used = (uint8_t *)calloc((size_t)tp->cache.size, 1);
}

static void tp_free(TotientPrefix *tp) {
    free(tp->pref0);
    free(tp->pref1);
    free(tp->pref2);
    free(tp->cache.entries);
    free(tp->cache.used);
}

/* ---- Problem-specific functions ---- */
static i64 nonconcurrent_candidate_count(i64 m, i64 n) {
    i64 mm = m % MOD;
    i64 mm1 = (mm - 1 + MOD) % MOD;
    i64 nn = n % MOD;
    i64 nn1 = (nn - 1 + MOD) % MOD;
    i64 nn1p = (nn + 1) % MOD;

    i64 two_same = mulmod(mulmod(mulmod(mulmod(mulmod(mm, mm1), nn), nn1), nn1p), INV6);

    i64 diff = (c3_mod(n + 2) - (n % MOD) + MOD) % MOD;
    i64 distinct = mulmod(c3_mod(m), diff);

    return (two_same + distinct) % MOD;
}

static i64 weighted_gcd_sum(i64 m, i64 n, TotientPrefix *tp) {
    i64 m1 = m - 1;
    i64 n1 = n - 1;
    i64 upper = m1 < n1 ? m1 : n1;
    i64 total = 0;

    i64 l = 1;
    while (l <= upper) {
        i64 qm = m1 / l;
        i64 qn = n1 / l;
        i64 r = m1 / qm;
        if (n1 / qn < r) r = n1 / qn;
        if (upper < r) r = upper;

        i64 r0, r1, r2, l0, l1, l2;
        tp_values(tp, r, &r0, &r1, &r2);
        tp_values(tp, l - 1, &l0, &l1, &l2);

        i64 s0 = (r0 - l0 + MOD) % MOD;
        i64 s1 = (r1 - l1 + MOD) % MOD;
        i64 s2 = (r2 - l2 + MOD) % MOD;

        i64 qm_mod = qm % MOD;
        i64 qn_mod = qn % MOD;

        i64 a0m = mulmod(qm_mod, m % MOD);
        i64 a1m = mulmod(mulmod(-qm_mod, (qm + 1) % MOD), INV2);
        i64 a0n = mulmod(qn_mod, n % MOD);
        i64 a1n = mulmod(mulmod(-qn_mod, (qn + 1) % MOD), INV2);

        i64 c0 = mulmod(a0m, a0n);
        i64 c1 = (mulmod(a0m, a1n) + mulmod(a1m, a0n)) % MOD;
        i64 c2 = mulmod(a1m, a1n);

        i128 sum = (i128)c0 * s0 + (i128)c1 * s1 + (i128)c2 * s2;
        i64 sum_mod = (i64)(sum % MOD);
        if (sum_mod < 0) sum_mod += MOD;
        total = (total + sum_mod) % MOD;

        l = r + 1;
    }
    return total;
}

static i64 concurrent_triple_count(i64 m, i64 n, TotientPrefix *tp) {
    i64 gcd_part = weighted_gcd_sum(m, n, tp);
    i64 endpoint = mulmod(c2_mod(m), c2_mod(n));
    return (gcd_part - endpoint + MOD) % MOD;
}

long long p994_native(void) {
    INV2 = (MOD + 1) / 2;
    INV6 = mod_pow(6, MOD - 2);

    i64 sieve_limit = 10000000;
    TotientPrefix tp;
    tp_init(&tp, sieve_limit);

    i64 m = 1234LL * 100000000LL;
    i64 n = 2345LL * 100000000LL;

    i64 nc = nonconcurrent_candidate_count(m, n);
    i64 ct = concurrent_triple_count(m, n, &tp);
    i64 result = (nc - ct + MOD) % MOD;

    tp_free(&tp);
    return result;
}
