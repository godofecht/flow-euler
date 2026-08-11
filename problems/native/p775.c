// Project Euler 775: Saving Paper
//
// Model n unit cubes as a polycube. With full-contact wrapping the paper
// equals the surface area. s(n) is the minimum surface area over all
// polycubes of n cubes; g(n) = 6n - s(n) is the paper saved. We need
// G(N) = sum_{n=1..N} g(n) mod 1e9+7 with N = 10^16.
//
// Surface-minimizing polycubes have a closed structure: between k^3 and
// (k+1)^3 cubes one starts from a k*k*k block and adds three orthogonal
// layers of sizes k^2, k(k+1), (k+1)^2. Within a layer only O(sqrt(m)) of
// the first m added cubes raise the surface by 2; the rest raise it by 0
// (the first cube of each layer raises it by 4). This gives a closed form
// for s(n) and for summing s(n) over huge ranges.
//
// Ported from the reference Python solver.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define MOD 1000000007LL

/* floor(cuberoot(x)) for x >= 0 via integer binary search. */
static long long icbrt_floor(unsigned long long x) {
    if (x == 0) return 0;
    unsigned long long lo = 0, hi = 1;
    while (hi * hi * hi <= x) hi *= 2;
    while (lo + 1 < hi) {
        unsigned long long mid = lo + (hi - lo) / 2;
        if (mid * mid * mid <= x) lo = mid;
        else hi = mid;
    }
    return (long long)lo;
}

/* floor(sqrt(x)) for x >= 0. */
static long long isqrt_ll(unsigned long long x) {
    if (x == 0) return 0;
    if (x < 2) return 1;
    unsigned long long r = x, t = x / 2;
    while (t < r) { r = t; t = (r + x / r) / 2; }
    return (long long)r;
}

/* Count of 'turn' positions in the first m steps of the square spiral layer:
   c(m) = max(0, isqrt(4m-1) - 1). */
static long long c_count(long long m) {
    if (m <= 1) return 0;
    return isqrt_ll((unsigned long long)(4 * m - 1)) - 1;
}

/* F(t) = sum_{m=1..t} c_count(m), O(1) via block structure of floor(sqrt). */
static long long c_prefix_sum(long long t) {
    if (t <= 1) return 0;
    long long a = isqrt_ll((unsigned long long)t); /* a^2 <= t < (a+1)^2 */

    /* sum_{i=1..a-1} (4i^2 + i) = (a-1)*a*(8a-1)/6 */
    long long base = (a - 1) * a * (8 * a - 1) / 6;

    long long sq = a * a;
    if (t == sq) return base;

    /* m in [a^2+1, a^2+a]: c(m) = 2a-1 */
    long long end_even = t < sq + a ? t : sq + a;
    long long cnt_even = end_even - (sq + 1) + 1;
    long long partial = cnt_even * (2 * a - 1);

    /* m in [a^2+a+1, (a+1)^2]: c(m) = 2a */
    if (t > sq + a) {
        long long cnt_odd = t - (sq + a + 1) + 1;
        partial += cnt_odd * (2 * a);
    }
    return base + partial;
}

static long long mod128(long long a, long long b, long long m) {
    /* (a*b) mod m via __int128 */
    return (long long)((__int128)a * (__int128)b % m);
}

/* S(N) = sum_{n=1..N} smin(n) mod mod. */
static long long sum_smin_mod(long long N, long long mod) {
    if (N <= 0) return 0;
    long long total = 6 % mod; /* n=1 */
    if (N == 1) return total;

    long long k_max = icbrt_floor((unsigned long long)(N - 1));
    for (long long k = 1; k <= k_max; k++) {
        long long k3 = k * k * k;
        if (k3 + 1 > N) break;
        long long full = (k + 1) * (k + 1) * (k + 1) - k3;
        long long L = N - k3;
        if (L > full) L = full;

        long long k2 = k * k;
        long long cap2 = k * (k + 1);

        long long lenA = L < k2 ? L : k2;
        long long rem = L - lenA;
        long long lenB = 0, lenC = 0;
        if (rem > 0) {
            lenB = rem < cap2 ? rem : cap2;
            rem -= lenB;
            if (rem > 0) lenC = rem;
        }

        long long c_k2 = c_count(k2);
        long long c_kk1 = c_count(cap2);

        /* sum of (c(pz)+c(qz)+c(rz)) across the block */
        long long sum_c = c_prefix_sum(lenA);
        if (lenB) sum_c += lenB * c_k2 + c_prefix_sum(lenB);
        if (lenC) sum_c += lenC * (c_k2 + c_kk1) + c_prefix_sum(lenC);

        /* sum of per-n constant term bv (4/8/12) across the block */
        long long sum_bv = 4 * lenA + 8 * lenB + 12 * lenC;

        /* block = L*(6*k2) + sum_bv + 2*sum_c, all mod.
           Use __int128 to avoid overflow on the big product. */
        long long term1 = mod128(L % mod, (6 * k2) % mod, mod);
        long long term2 = (sum_bv % mod + mod128(2 % mod, sum_c % mod, mod)) % mod;
        long long block = (term1 + term2) % mod;
        total = (total + block) % mod;
    }
    return total;
}

/* G(N) = sum_{n=1..N} g(n) mod mod, g(n) = 6n - smin(n). */
static long long G_mod(long long N, long long mod) {
    /* sum_{n=1..N} 6n = 3N(N+1) */
    long long sum6 = mod128(mod128(3, N % mod, mod), ((N + 1) % mod), mod);
    long long S = sum_smin_mod(N, mod);
    long long r = (sum6 - S) % mod;
    if (r < 0) r += mod;
    return r;
}

long long p775_native(void) {
    /* Self-test against known values. */
    /* g(10) = 30, g(18) = 66, G(18) = 530, G(10^6) = 951640919 */
    return G_mod(10000000000000000LL, MOD); /* N = 10^16 */
}
