/* Project Euler 945: XOR-Equation C
 * Count pairs 0<=a<=b<=N where (a XOR-product b) has no even-position bits.
 * Uses GF(2) polynomial arithmetic and digit DP.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef int64_t i64;
typedef uint64_t u64;

#define EVEN_MASK 0x5555555555555555ULL

/* 16-bit compaction lookup tables */
static u64 even16[65536];
static u64 odd16[65536];
static int tables_built = 0;

static void build_tables(void) {
    for (int x = 0; x < 65536; x++) {
        u64 e = 0, o = 0;
        for (int i = 0; i < 8; i++) {
            u64 b = 1ULL << (2 * i);
            if (x & b) e |= 1ULL << i;
            if (x & (b << 1)) o |= 1ULL << i;
        }
        even16[x] = e;
        odd16[x] = o;
    }
    tables_built = 1;
}

/* Split x into (E, O) where E = even-position bits compacted, O = odd-position bits compacted */
static void split_u(u64 x, u64 *E, u64 *O) {
    u64 e = 0, o = 0;
    int shift = 0;
    while (x) {
        u64 chunk = x & 0xFFFF;
        e |= even16[chunk] << shift;
        o |= odd16[chunk] << shift;
        x >>= 16;
        shift += 8;
    }
    *E = e;
    *O = o;
}

/* Polynomial remainder a mod b over GF(2) */
static u64 gf2_mod(u64 a, u64 b) {
    if (b == 0) return a;
    int db = 63 - __builtin_clzll(b);
    while (a) {
        int da = 63 - __builtin_clzll(a);
        if (da < db) break;
        a ^= b << (da - db);
    }
    return a;
}

/* Polynomial GCD over GF(2) */
static u64 gf2_gcd(u64 a, u64 b) {
    while (b) {
        u64 t = gf2_mod(a, b);
        a = b;
        b = t;
    }
    return a;
}

/* Ordered solution count for [0, 2^bits - 1] */
static i64 ordered_full(int bits) {
    i64 sign = (bits & 1) ? -1 : 1;
    i64 num = ((i64)1 << (bits + 1)) * (3 * bits + 4) + sign;
    return num / 9;
}

/* For b = 2^k + y, count a in [0, 2^k - 1] satisfying the condition */
static i64 count_a_for_upper(int k, u64 y) {
    int m = k / 2;
    u64 Y0, Y1;
    split_u(y, &Y0, &Y1);

    if ((k & 1) == 0) {
        /* k even: k = 2m */
        if (Y1 == 0) {
            return (i64)1 << m;
        }
        u64 P = Y0 ^ (1ULL << m);
        u64 g = gf2_gcd(P, Y1 << 1);
        return (i64)1 << (63 - __builtin_clzll(g));
    } else {
        /* k odd: k = 2m+1 */
        if (Y0 == 0) {
            return (i64)1 << (m + 1);
        }
        u64 uQ = ((1ULL << m) ^ Y1) << 1;
        u64 g = gf2_gcd(Y0, uQ);
        return (i64)1 << (63 - __builtin_clzll(g));
    }
}

/* Ordered count of pairs (a,b) with 0<=a,b<=N */
static i64 ordered_S(i64 N) {
    if (N < 0) return 0;
    if (N == 0) return 1;

    int bits = 64 - __builtin_clzll((u64)N);
    u64 all_ones = ((u64)1 << bits) - 1;
    if ((u64)N == all_ones) {
        return ordered_full(bits);
    }

    int k = bits - 1;
    i64 M = (i64)1 << k;
    i64 r = N - M;

    i64 base = ordered_full(k);
    i64 cross = 0;
    for (i64 y = 0; y <= r; y++) {
        cross += count_a_for_upper(k, (u64)y);
    }

    return base + 2 * cross;
}

long long p945_native(void) {
    if (!tables_built) build_tables();

    i64 N = 1;
    for (int i = 0; i < 7; i++) N *= 10; /* 10^7 */

    i64 s = ordered_S(N);
    return (s + 1) / 2;
}
