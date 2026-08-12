// Project Euler 889: Rational Blancmange.
// F(k, t, r) = (2^(2k) - 1) * T(((2^t + 1)^r) / (2^k + 1)) mod 1000062031.
// k = 10^18 + 31, t = 10^14 + 31, r = 62.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned long long u64;
typedef long long i64;
typedef __int128 i128;

#define MOD 1000062031ULL

static u64 mulmod(u64 a, u64 b, u64 m) {
    return (u64)((i128)a * b % m);
}

static u64 powmod(u64 a, u64 e, u64 m) {
    u64 r = 1 % m; a %= m;
    while (e > 0) {
        if (e & 1) r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}

// Count set bits
static int popcount_u64(u64 x) { return __builtin_popcountll(x); }

// Bit length
static int bitlen_u64(u64 x) {
    if (x == 0) return 0;
    return 64 - __builtin_clzll(x);
}

// Compute binomial C(n, k) for small n (up to 62)
static u64 binom_small(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    u64 r = 1;
    for (int i = 0; i < k; i++) {
        r = r * (u64)(n - i) / (u64)(i + 1);
    }
    return r;
}

// Brute F_mod for small k (Q = 2^k + 1 fits in u64, so k <= 63)
static u64 brute_F_mod(u64 k, u64 t, u64 r, u64 mod) {
    u64 Q = (1ULL << k) + 1;
    u64 base = (1ULL << t) + 1;
    u64 N_mod = powmod(base, r, Q);
    u64 m = N_mod;
    u64 inv2 = powmod(2, mod - 2, mod);
    u64 weight = powmod(2, k, mod);
    u64 ans = 0;
    for (u64 j = 0; j < k; j++) {
        u64 d = (m <= Q - m) ? m : (Q - m);
        ans = (ans + mulmod(d % mod, weight, mod)) % mod;
        weight = mulmod(weight, inv2, mod);
        m = mulmod(m, 2, Q);
    }
    return ans;
}

// Compute bit positions of N = (2^t + 1)^r using sparse binomial expansion.
// Returns number of positions, fills positions array (must be pre-allocated).
static int bit_positions_of_N(u64 t, int r, u64 *positions) {
    // max_coeff_bitlen = bit_length(C(r, r/2))
    int half = r / 2;
    u64 max_coeff = binom_small(r, half);
    int max_coeff_bitlen = bitlen_u64(max_coeff);

    int n = 0;
    for (int u = 0; u <= r; u++) {
        u64 c = binom_small(r, u);
        u64 base = t * (u64)u;
        // For each set bit b in c, add position base + b
        while (c > 0) {
            int b = bitlen_u64(c & (~(c - 1))) - 1;  // position of lowest set bit
            // Actually, let's use __builtin_ctzll
            b = __builtin_ctzll(c);
            positions[n++] = base + (u64)b;
            c &= c - 1;  // clear lowest set bit
        }
    }
    // Sort positions (simple insertion sort since n is small, ~3800 max)
    for (int i = 1; i < n; i++) {
        u64 key = positions[i];
        int j = i - 1;
        while (j >= 0 && positions[j] > key) {
            positions[j + 1] = positions[j];
            j--;
        }
        positions[j + 1] = key;
    }
    return n;
}

// Fast F_mod for large k
static u64 fast_F_mod(u64 k, u64 t, int r, u64 mod) {
    u64 positions[4000];
    int n = bit_positions_of_N(t, r, positions);
    u64 max_pos = positions[n - 1];

    u64 pow2_k = powmod(2, k, mod);

    // Precompute 2^p mod mod and 2^{k+p} mod mod for each position
    u64 *vals_low = (u64 *)malloc(n * sizeof(u64));
    u64 *vals_high = (u64 *)malloc(n * sizeof(u64));
    for (int i = 0; i < n; i++) {
        vals_low[i] = powmod(2, positions[i], mod);
        vals_high[i] = mulmod(pow2_k, vals_low[i], mod);
    }

    // Prefix sums
    u64 *prefix_low = (u64 *)malloc((n + 1) * sizeof(u64));
    u64 *prefix_high = (u64 *)malloc((n + 1) * sizeof(u64));
    prefix_low[0] = 0;
    prefix_high[0] = 0;
    for (int i = 0; i < n; i++) {
        prefix_low[i + 1] = (prefix_low[i] + vals_low[i]) % mod;
        prefix_high[i + 1] = (prefix_high[i] + vals_high[i]) % mod;
    }
    u64 total_low = prefix_low[n];

    // Base sum: for each bit p, contribution = (k-p)*2^{k+p} - p*2^p
    u64 ans = 0;
    for (int i = 0; i < n; i++) {
        u64 p = positions[i];
        u64 v_low = vals_low[i];
        u64 v_high = vals_high[i];
        // (k - p) * v_high - p * v_low
        u64 km = k % mod;
        u64 pm = p % mod;
        u64 term = mulmod((km + mod - pm) % mod, v_high, mod);
        term = (term + mod - mulmod(pm, v_low, mod)) % mod;
        ans = (ans + term) % mod;
    }

    // Corrections: for i = 1..n-1, p0 = positions[i]
    for (int i = 1; i < n; i++) {
        u64 p0 = positions[i];
        // S = sum_high_le - sum_low_gt
        u64 sum_high_le = prefix_high[i + 1];
        u64 sum_low_gt = (total_low + mod - prefix_low[i + 1]) % mod;
        u64 S = (sum_high_le + mod - sum_low_gt) % mod;

        // pow2_p0_plus1 = 2^{p0+1} mod mod = vals_low[i] * 2
        u64 pow2_p0_plus1 = mulmod(vals_low[i], 2, mod);
        // q_times = (2^k + 1) * 2^{p0+1} mod mod
        u64 q_times = mulmod((pow2_k + 1) % mod, pow2_p0_plus1, mod);
        // delta = q_times - 2*S
        u64 delta = (q_times + mod - mulmod(S, 2, mod)) % mod;
        ans = (ans + delta) % mod;
    }

    free(vals_low);
    free(vals_high);
    free(prefix_low);
    free(prefix_high);
    return ans;
}

static u64 F_mod(u64 k, u64 t, u64 r, u64 mod) {
    if (k <= 63) {
        return brute_F_mod(k, t, r, mod);
    }
    return fast_F_mod(k, t, (int)r, mod);
}

long long p889_native(void) {
    // Test cases
    // F_mod(3, 1, 1, MOD) == 42 -- uses brute (k=3 <= 63)
    // F_mod(13, 3, 3, MOD) == 23093880 -- uses brute (k=13 <= 63)
    // F_mod(103, 13, 6, MOD) == 878922518 -- uses fast (k=103 > 63)

    u64 k = 1000000000000000031ULL;  // 10^18 + 31
    u64 t = 100000000000031ULL;      // 10^14 + 31
    u64 r = 62;

    return (long long)F_mod(k, t, r, MOD);
}
