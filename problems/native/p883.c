// Project Euler 883: Remarkable Triangles on the triangular/hexagonal lattice.
// Counts triangles with lattice-point incenter and inradius <= R.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef unsigned long long u64;
typedef long long i64;
typedef uint32_t u32;
typedef uint8_t u8;
typedef __int128 i128;

static u64 isqrt_u64(u64 n) {
    if (n == 0) return 0;
    u64 x = (u64)sqrtl((long double)n);
    while (x * x > n) x--;
    while ((x + 1) * (x + 1) <= n) x++;
    return x;
}

static u64 gcd_u64(u64 a, u64 b) {
    while (b) { u64 t = a % b; a = b; b = t; }
    return a;
}

// ---- Hash table for cache: B -> count ----
#define CACHE_CAP (1 << 20)  // 1M slots
#define CACHE_MASK (CACHE_CAP - 1)

typedef struct {
    u64 key;
    u64 val;
    u8 used;
} CacheEntry;

static CacheEntry g_cache[CACHE_CAP];

static u64 cache_hash(u64 key) {
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return key & CACHE_MASK;
}

static u64 cache_lookup(u64 key, int *found) {
    u64 h = cache_hash(key);
    while (g_cache[h].used) {
        if (g_cache[h].key == key) { *found = 1; return g_cache[h].val; }
        h = (h + 1) & CACHE_MASK;
    }
    *found = 0;
    return 0;
}

static void cache_insert(u64 key, u64 val) {
    u64 h = cache_hash(key);
    while (g_cache[h].used) {
        if (g_cache[h].key == key) { g_cache[h].val = val; return; }
        h = (h + 1) & CACHE_MASK;
    }
    g_cache[h].used = 1;
    g_cache[h].key = key;
    g_cache[h].val = val;
}

// prefix sum of chi on 1..m: count of 1 mod 3 minus count of 2 mod 3
static i64 chi_prefix(i64 m) {
    return (m + 2) / 3 - (m + 1) / 3;
}

static u64 count_hex_points_leq(u64 B) {
    if (B == 0) return 1;
    int found;
    u64 hit = cache_lookup(B, &found);
    if (found) return hit;

    i64 total = 0;
    i64 n = (i64)B;
    i64 i = 1;
    while (i <= n) {
        i64 q = n / i;
        i64 j = n / q;
        total += q * (chi_prefix(j) - chi_prefix(i - 1));
        i = j + 1;
    }

    u64 res = 1 + 6 * (u64)total;
    cache_insert(B, res);
    return res;
}

static u64 remarkable_triangles(u64 R_num, u64 R_den) {
    u64 C_num = 12 * R_num * R_num;
    u64 C_den = R_den * R_den;
    u64 Dmax = isqrt_u64(C_num / C_den);

    u32 *M = (u32 *)calloc(Dmax + 1, sizeof(u32));

    // Family 1: d = 3*u*v, gcd(u,v)=1, u>v, u mod 3 != v mod 3
    // Multiplicity = 2^{omega(t)-1} where t = u*v, t % 3 != 1
    u64 Nmax = Dmax / 3;
    u8 *omega = (u8 *)calloc(Nmax + 1, sizeof(u8));
    for (u64 p = 2; p <= Nmax; p++) {
        if (omega[p] == 0) {  // prime
            for (u64 k = p; k <= Nmax; k += p) {
                omega[k]++;
            }
        }
    }
    for (u64 t = 2; t <= Nmax; t++) {
        if (t % 3 == 1) continue;
        M[3 * t] += 1u << (omega[t] - 1);
    }
    free(omega);

    // Family 2: d = (u-v)(u+2v) = a*(a+3v), gcd(u,v)=1, u mod 3 != v mod 3
    u64 D = Dmax;
    u64 vmax = D / 3;
    for (u64 v = 1; v <= vmax; v++) {
        // a^2 + 3v*a - D <= 0 => a <= floor((sqrt(9v^2+4D) - 3v)/2)
        u64 disc = 9 * v * v + 4 * D;
        u64 sq = isqrt_u64(disc);
        i64 amax = ((i64)sq - (i64)(3 * v)) / 2;
        if (amax <= 0) continue;
        for (i64 a = 1; a <= amax; a++) {
            u64 u = v + (u64)a;
            if (u % 3 == v % 3) continue;
            if (gcd_u64(u, v) != 1) continue;
            u64 d = (u64)a * ((u64)a + 3 * v);
            M[d]++;
        }
    }

    // Sum contributions over d
    u64 total_scalene = 0;
    for (u64 d = 1; d <= Dmax; d++) {
        u32 mult = M[d];
        if (mult == 0) continue;

        u64 B = C_num / (C_den * d * d);
        if (B == 0) continue;

        u64 pts;
        if (d % 3 == 0) {
            pts = count_hex_points_leq(B) - 1;
        } else {
            pts = count_hex_points_leq(B / 3) - 1;
        }
        total_scalene += 2 * (u64)mult * pts;
    }

    // Equilateral triangles
    u64 Beq = (4 * R_num * R_num) / (R_den * R_den);
    u64 equi_points = count_hex_points_leq(Beq) - 1;
    u64 total_equilateral = equi_points / 3;

    free(M);
    return total_scalene + total_equilateral;
}

long long p883_native(void) {
    return (long long)remarkable_triangles(1000000, 1);
}
