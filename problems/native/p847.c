// Project Euler 847: Jack's Bean
// H(R_19) mod 1e9+7, where R_19 = repunit(19)
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

#define MOD 1000000007LL

static i64 modpow(i64 base, i64 exp, i64 mod) {
    i64 r = 1; base %= mod; if (base < 0) base += mod;
    while (exp > 0) {
        if (exp & 1) r = (i64)((i128)r * base % mod);
        base = (i64)((i128)base * base % mod);
        exp >>= 1;
    }
    return r;
}

static i64 INV2, INV6;

// ceil_log2(n) for n >= 1, 0 for n <= 1
static i64 ceil_log2(i64 n) {
    if (n <= 1) return 0;
    return 64 - __builtin_clzll((u64)(n - 1));
}

// repunit(n) = 111...1 (n digits)
static i64 repunit(int n) {
    i64 r = 0;
    for (int i = 0; i < n; i++) r = r * 10 + 1;
    return r;
}

// ---- Modular sum helpers ----

static i64 sum_1(i64 n) {
    if (n <= 0) return 0;
    i64 a = n % MOD, b = (n + 1) % MOD;
    return (i64)((i128)a * b % MOD * INV2 % MOD);
}

static i64 sum_2(i64 n) {
    if (n <= 0) return 0;
    i64 a = n % MOD, b = (n + 1) % MOD, c = (2 * n + 1) % MOD;
    return (i64)((i128)((i128)a * b % MOD) * c % MOD * INV6 % MOD);
}

static i64 range_sum_1(i64 l, i64 r) {
    if (l > r) return 0;
    i64 s = sum_1(r) - sum_1(l - 1);
    return (s % MOD + MOD) % MOD;
}

static i64 range_sum_2(i64 l, i64 r) {
    if (l > r) return 0;
    i64 s = sum_2(r) - sum_2(l - 1);
    return (s % MOD + MOD) % MOD;
}

static i64 sum_triples_count(i64 l, i64 r) {
    if (l > r) return 0;
    i64 cnt = (r - l + 1) % MOD;
    i64 s1 = range_sum_1(l, r);
    i64 s2 = range_sum_2(l, r);
    i64 total = (s2 + 3 * s1 % MOD + 2 * cnt % MOD) % MOD;
    return (i64)((i128)total * INV2 % MOD);
}

// ---- Bad count helpers ----

static i64 threshold_t(int k) {
    if (k < 3) return (i64)1e18; // effectively infinity
    return 3LL * (1LL << (k - 2)) + 2; // 3 * 2^(k-2) + 2
}

static i64 base_bad_for_block(i64 M, i64 L) {
    if (L < (M >> 1) + 2) return 0;
    i64 n = 2 * L - M - 1;
    return n * (n - 1) / 2;
}

// Precomputed bad_count_for_sum for s = 1..16
static i64 bad_count_small[17];

static i64 bad_count_for_sum(i64 s) {
    if (s <= 7) return 0;
    // Python: k = (s - 1).bit_length() - 1, s in (2^k, 2^(k+1)]
    int k = 64 - __builtin_clzll((u64)(s - 1)) - 1;
    i64 M = 1LL << k;
    i64 L = s - M;
    i64 res = base_bad_for_block(M, L);
    if (k >= 3 && L >= threshold_t(k)) {
        res += 3 * bad_count_for_sum(L);
    }
    return res;
}

static void init_bad_count_small(void) {
    for (int s = 0; s <= 16; s++)
        bad_count_small[s] = bad_count_for_sum(s);
}

// ---- sum_base_bad_block_mod ----

static i64 sum_base_bad_block_mod(i64 M, i64 R) {
    if (R <= 0) return 0;
    i64 start = (M >> 1) + 2;
    if (R < start) return 0;
    i64 a = start, b = R;
    i64 cnt = (b - a + 1) % MOD;
    i64 sumL = range_sum_1(a, b);
    i64 sumL2 = range_sum_2(a, b);
    i64 const_val = (i64)((i128)((M + 1) % MOD) * ((M + 2) % MOD) % MOD * INV2 % MOD);
    i64 res = (2 * sumL2 % MOD - (i64)((i128)((2 * M + 3) % MOD) * sumL % MOD) + (i64)((i128)const_val * cnt % MOD)) % MOD;
    return (res % MOD + MOD) % MOD;
}

// ---- Memoization hash map for bad_prefix_sum_mod ----

#define HT_CAP 4096

typedef struct {
    i64 key;
    i64 val;
    char used;
} HTEntry;

static HTEntry ht[HT_CAP];

static void ht_clear(void) {
    memset(ht, 0, sizeof(ht));
}

static u64 ht_hash(i64 key) {
    return (u64)key * 1099511628211ULL;
}

static i64 *ht_lookup(i64 key) {
    u64 h = ht_hash(key) % HT_CAP;
    while (ht[h].used) {
        if (ht[h].key == key) return &ht[h].val;
        h = (h + 1) % HT_CAP;
    }
    return NULL;
}

static void ht_insert(i64 key, i64 val) {
    u64 h = ht_hash(key) % HT_CAP;
    while (ht[h].used) {
        if (ht[h].key == key) { ht[h].val = val; return; }
        h = (h + 1) % HT_CAP;
    }
    ht[h].used = 1;
    ht[h].key = key;
    ht[h].val = val;
}

// Forward declarations
static i64 bad_prefix_sum_mod(i64 N);
static i64 block_bad_sum_mod(i64 M, i64 R);

static i64 block_bad_sum_mod(i64 M, i64 R) {
    if (R <= 0) return 0;
    int k = 63 - __builtin_clzll((u64)M); // log2(M)
    i64 res = sum_base_bad_block_mod(M, R);
    i64 T = threshold_t(k);
    if (T <= R) {
        i64 diff = (bad_prefix_sum_mod(R) - bad_prefix_sum_mod(T - 1)) % MOD;
        res = (res + 3 * diff % MOD) % MOD;
    }
    return (res % MOD + MOD) % MOD;
}

static i64 bad_prefix_sum_mod(i64 N) {
    if (N <= 0) return 0;
    if (N <= 16) {
        i64 total = 0;
        for (int s = 1; s <= N; s++) total += bad_count_small[s];
        return total % MOD;
    }

    // Check memoization
    i64 *cached = ht_lookup(N);
    if (cached) return *cached;

    i64 pow2 = 1LL << (63 - __builtin_clzll((u64)N));
    i64 result;
    if (N == pow2) {
        if (pow2 == 1) {
            result = 0;
        } else {
            i64 M = pow2 >> 1;
            result = (bad_prefix_sum_mod(M) + block_bad_sum_mod(M, M)) % MOD;
        }
    } else {
        result = (bad_prefix_sum_mod(pow2) + block_bad_sum_mod(pow2, N - pow2)) % MOD;
    }

    ht_insert(N, result);
    return result;
}

// ---- base_part_mod ----

static i64 base_part_mod(i64 N) {
    if (N <= 0) return 0;
    i64 res = 0;
    i64 m = 1;
    i64 low = 2;
    while (low <= N) {
        i64 high = (1LL << m);
        if (high > N) high = N;
        res = (res + (i64)((i128)m * sum_triples_count(low, high) % MOD)) % MOD;
        m++;
        low = (1LL << (m - 1)) + 1;
    }
    return res;
}

// ---- H_mod ----

static i64 H_mod(i64 N) {
    return (base_part_mod(N) + bad_prefix_sum_mod(N)) % MOD;
}

long long p847_native(void) {
    INV2 = (MOD + 1) / 2;
    INV6 = modpow(6, MOD - 2, MOD);

    init_bad_count_small();
    ht_clear();

    i64 n = repunit(19);
    return H_mod(n);
}
