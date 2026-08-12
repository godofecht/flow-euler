// Project Euler 858: LCM
// G(800) mod 1e9+7, where G(N) = sum over all subsets S of {1..N} of lcm(S)
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

#define MOD 1000000007LL
#define N 800
#define SQRT_N 28
#define N_WORDS ((N + 63) / 64)  // 13
#define KMAX (N / (SQRT_N + 1))  // 800/29 = 27

static i64 modpow(i64 base, i64 exp, i64 mod) {
    i64 r = 1; base %= mod; if (base < 0) base += mod;
    while (exp > 0) {
        if (exp & 1) r = (i64)((i128)r * base % mod);
        base = (i64)((i128)base * base % mod);
        exp >>= 1;
    }
    return r;
}

// Bitset for maskN (N bits)
typedef struct { uint64_t bits[N_WORDS]; } Bitset;

static void bitset_or(Bitset *dst, const Bitset *src) {
    for (int i = 0; i < N_WORDS; i++) dst->bits[i] |= src->bits[i];
}

static int bitset_popcount(const Bitset *bs) {
    int cnt = 0;
    for (int i = 0; i < N_WORDS; i++) cnt += __builtin_popcountll(bs->bits[i]);
    return cnt;
}

static void bitset_clear(Bitset *bs) {
    memset(bs->bits, 0, sizeof(bs->bits));
}

// Sieve
static int is_prime[N + 1];
static int primes[N + 1];
static int num_primes;

static void sieve(void) {
    memset(is_prime, 1, sizeof(is_prime));
    is_prime[0] = is_prime[1] = 0;
    for (int p = 2; p * p <= N; p++) {
        if (is_prime[p]) {
            for (int m = p * p; m <= N; m += p) is_prime[m] = 0;
        }
    }
    num_primes = 0;
    for (int i = 2; i <= N; i++) if (is_prime[i]) primes[num_primes++] = i;
}

// Global state
static i64 e_val[N + 1];        // e[p] = max exponent
static i64 p_pow_e_mod[N + 1];  // p^e mod MOD
static i64 inv_p_pow_e[N + 1];  // inverse of p^e mod MOD
static i64 pow2[N + 1];         // 2^i mod MOD
static i64 inv2pow[KMAX + 1];   // (1/2)^i mod MOD
static i64 prefix_masks[KMAX + 1]; // (1 << t) - 1

static int small_primes[32];
static int num_small_primes;
static int large_primes[200];
static int num_large_primes;

static i64 w_large[N + 1];  // (p-1)/p mod MOD for large primes

// Options for each small prime
typedef struct {
    Bitset maskN;
    uint32_t maskK;
    i64 weight;
} Option;

static Option options[32][16];  // [small_prime_index][option_index]
static int num_options[32];     // number of options per small prime

// Hash map for large_product memoization
#define LP_CAP 65536
typedef struct {
    uint32_t key;
    i64 val;
    char used;
} LPEntry;

static LPEntry lp_table[LP_CAP];

static void lp_clear(void) { memset(lp_table, 0, sizeof(lp_table)); }

static i64 *lp_lookup(uint32_t key) {
    u64 h = (u64)key * 2654435761ULL % LP_CAP;
    while (lp_table[h].used) {
        if (lp_table[h].key == key) return &lp_table[h].val;
        h = (h + 1) % LP_CAP;
    }
    return NULL;
}

static void lp_insert(uint32_t key, i64 val) {
    u64 h = (u64)key * 2654435761ULL % LP_CAP;
    while (lp_table[h].used) {
        if (lp_table[h].key == key) { lp_table[h].val = val; return; }
        h = (h + 1) % LP_CAP;
    }
    lp_table[h].used = 1;
    lp_table[h].key = key;
    lp_table[h].val = val;
}

static i64 large_product(uint32_t maskK) {
    i64 *cached = lp_lookup(maskK);
    if (cached) return *cached;

    i64 prod = 1;
    for (int i = 0; i < num_large_primes; i++) {
        int p = large_primes[i];
        int t = N / p;  // <= KMAX
        int covered = 0;
        if (t > 0) {
            covered = __builtin_popcount(maskK & (uint32_t)prefix_masks[t]);
        }
        int new_val = t - covered;
        // factor = 1 - w_large[p] * inv2pow[new]
        i64 factor = (1 - (i64)((i128)w_large[p] * inv2pow[new_val] % MOD)) % MOD;
        factor = (factor % MOD + MOD) % MOD;
        prod = (i64)((i128)prod * factor % MOD);
    }
    lp_insert(maskK, prod);
    return prod;
}

static i64 dfs_total;

static void dfs(int i, Bitset maskN, uint32_t maskK, i64 coeff) {
    if (i == num_small_primes) {
        int covered_count = bitset_popcount(&maskN);
        i64 base = pow2[N - covered_count];
        i64 lp = large_product(maskK);
        dfs_total = (dfs_total + (i64)((i128)coeff * base % MOD * lp % MOD)) % MOD;
        return;
    }
    for (int j = 0; j < num_options[i]; j++) {
        Bitset new_maskN = maskN;
        bitset_or(&new_maskN, &options[i][j].maskN);
        uint32_t new_maskK = maskK | options[i][j].maskK;
        i64 new_coeff = (i64)((i128)coeff * options[i][j].weight % MOD);
        dfs(i + 1, new_maskN, new_maskK, new_coeff);
    }
}

long long p858_native(void) {
    sieve();

    // Compute e[p], p^e mod MOD, inv_p^e
    i64 L_mod = 1;
    for (int pi = 0; pi < num_primes; pi++) {
        int p = primes[pi];
        int ep = 0;
        i64 t = p;
        while (t <= N) { ep++; t *= p; }
        e_val[p] = ep;
        i64 pe_mod = modpow(p, ep, MOD);
        p_pow_e_mod[p] = pe_mod;
        inv_p_pow_e[p] = modpow(pe_mod, MOD - 2, MOD);
        L_mod = (i64)((i128)L_mod * pe_mod % MOD);
    }

    // Split primes
    num_small_primes = 0;
    num_large_primes = 0;
    for (int pi = 0; pi < num_primes; pi++) {
        int p = primes[pi];
        if (p <= SQRT_N) small_primes[num_small_primes++] = p;
        else large_primes[num_large_primes++] = p;
    }

    // Precompute pow2 and inv2pow
    pow2[0] = 1;
    for (int i = 1; i <= N; i++) pow2[i] = (i64)((i128)pow2[i-1] * 2 % MOD);
    i64 inv2 = (MOD + 1) / 2;
    inv2pow[0] = 1;
    for (int i = 1; i <= KMAX; i++) inv2pow[i] = (i64)((i128)inv2pow[i-1] * inv2 % MOD);

    // prefix_masks
    for (int t = 0; t <= KMAX; t++) prefix_masks[t] = (1LL << t) - 1;

    // Build options for each small prime
    for (int si = 0; si < num_small_primes; si++) {
        int p = small_primes[si];
        int ep = e_val[p];
        int nopt = 0;

        // r=0: not selected
        bitset_clear(&options[si][nopt].maskN);
        options[si][nopt].maskK = 0;
        options[si][nopt].weight = 1;
        nopt++;

        // r=1..ep
        for (int r = 1; r <= ep; r++) {
            i64 q = 1;
            for (int j = 0; j < r; j++) q *= p;

            Bitset maskN;
            bitset_clear(&maskN);
            for (i64 m = q; m <= N; m += q) {
                int pos = (int)(m - 1);
                maskN.bits[pos / 64] |= (1ULL << (pos % 64));
            }

            uint32_t maskK = 0;
            if (q <= KMAX) {
                for (i64 m = q; m <= KMAX; m += q) {
                    maskK |= (1U << (int)(m - 1));
                }
            }

            // phi(p^r) = p^(r-1) * (p-1)
            i64 phi = (i64)((i128)modpow(p, r - 1, MOD) * (p - 1) % MOD);
            i64 w = (i64)((i128)phi * inv_p_pow_e[p] % MOD);
            i64 weight = (MOD - w) % MOD;  // negative sign

            options[si][nopt].maskN = maskN;
            options[si][nopt].maskK = maskK;
            options[si][nopt].weight = weight;
            nopt++;
        }
        num_options[si] = nopt;
    }

    // Compute w_large for large primes
    for (int i = 0; i < num_large_primes; i++) {
        int p = large_primes[i];
        i64 invp = modpow(p, MOD - 2, MOD);
        w_large[p] = (i64)((i128)(p - 1) * invp % MOD);
    }

    // Run DFS
    lp_clear();
    dfs_total = 0;
    Bitset empty_mask;
    bitset_clear(&empty_mask);
    dfs(0, empty_mask, 0, 1);

    return (i64)((i128)L_mod * dfs_total % MOD);
}
