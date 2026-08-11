// Project Euler 720: sum of ranks of the "unpredictable" permutations.
// Port of the reference Python solver.  Builds the doubling arrays for
// size 2^24 then computes the Lehmer-code rank modulo 1e9+7.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define MOD 1000000007LL

static long long mulmod(long long a, long long b) {
    return (a % MOD) * (b % MOD) % MOD;
}

// Build vals[] and codes[] for a target size that is a power of two >= 4.
static void build_arrays(long target_size, uint32_t **vals_out, uint32_t **codes_out) {
    uint32_t *vals = malloc((size_t)target_size * sizeof(uint32_t));
    uint32_t *codes = malloc((size_t)target_size * sizeof(uint32_t));
    if (!vals || !codes) { fprintf(stderr, "oom\n"); exit(1); }

    vals[0] = 1; vals[1] = 3; vals[2] = 2; vals[3] = 4;
    codes[0] = 0; codes[1] = 1; codes[2] = 0; codes[3] = 0;
    long size = 4;

    while (size < target_size) {
        long m = size;
        long n = m << 1;
        uint32_t *nv = malloc((size_t)n * sizeof(uint32_t));
        uint32_t *nc = malloc((size_t)n * sizeof(uint32_t));
        if (!nv || !nc) { fprintf(stderr, "oom\n"); exit(1); }

        for (long i = 0; i < m - 1; ++i) {
            uint32_t v = vals[i];
            nv[i] = (v << 1) - 1;
            nc[i] = (v - 1) + codes[i];
        }
        nv[m - 1] = 2;
        nc[m - 1] = 0;
        uint32_t v_last = vals[m - 1];
        nv[m] = (v_last << 1) - 1;
        nc[m] = (uint32_t)(m - 2);
        for (long j = 1; j < m; ++j) {
            uint32_t v = vals[j];
            nv[m + j] = v << 1;
            nc[m + j] = codes[j];
        }

        free(vals);
        free(codes);
        vals = nv;
        codes = nc;
        size = n;
    }

    *vals_out = vals;
    *codes_out = codes;
}

static long long rank_from_prev(const uint32_t *vals, const uint32_t *codes, long m) {
    long long rank = 0;
    long long fact = 1;
    long long step = 1;

    for (long j = m - 1; j > 0; --j) {
        long long l = codes[j];
        rank = (rank + mulmod(l, fact)) % MOD;
        fact = mulmod(fact, step);
        step += 1;
    }
    {
        long long l = m - 2;
        rank = (rank + mulmod(l, fact)) % MOD;
        fact = mulmod(fact, step);
        step += 1;
    }
    fact = mulmod(fact, step);
    step += 1;

    for (long i = m - 2; i >= 0; --i) {
        uint32_t v = vals[i];
        long long l = (long long)(v - 1) + (long long)codes[i];
        rank = (rank + mulmod(l, fact)) % MOD;
        fact = mulmod(fact, step);
        step += 1;
    }
    return rank;
}

static long long solve_power(int k) {
    if (k == 0) return 1;
    if (k == 1) return 1;
    if (k == 2) return 3;
    long target = 1L << (k - 1);
    uint32_t *vals, *codes;
    build_arrays(target, &vals, &codes);
    long long rank0 = rank_from_prev(vals, codes, target);
    free(vals);
    free(codes);
    return (rank0 + 1) % MOD;
}

long long p720_native(void) {
    return solve_power(25);
}
