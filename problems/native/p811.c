/* Project Euler 811: Bitwise Recursion
   H(t, r) = A((2^t + 1)^r) mod 1000062031
   with t = 10^14 + 31, r = 62.
   Port of the Python reference solver. */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef __int128 i128;

#define MOD 1000062031LL

static int popcount_i128(i128 x) {
    int cnt = 0;
    while (x) {
        x &= x - 1;
        cnt++;
    }
    return cnt;
}

/* Max bit-length among binomial coefficients C(r, k) for k=0..r */
static int max_binom_bitlen(int r) {
    i128 c = 1;
    int mx = 1;
    for (int k = 0; k <= r; k++) {
        int bl = 0;
        i128 tmp = c;
        while (tmp) { bl++; tmp >>= 1; }
        if (bl > mx) mx = bl;
        if (k < r) {
            c = c * (r - k) / (k + 1);
        }
    }
    return mx;
}

/* Positions of 1-bits in (2^t + 1)^r using the binomial block method.
   Assumes t >= max_binom_bitlen(r) so blocks don't overlap.
   Returns positions in increasing order. */
static int one_positions_via_binom(long long t, int r, long long *pos_out) {
    int cnt = 0;
    i128 c = 1; /* C(r, 0) */
    for (int k = 0; k <= r; k++) {
        i128 x = c;
        while (x) {
            i128 lsb = x & (-x);
            /* bit position = log2(lsb) */
            int bit = 0;
            i128 tmp = lsb;
            while (tmp > 1) { tmp >>= 1; bit++; }
            pos_out[cnt++] = (long long)k * t + bit;
            x -= lsb;
        }
        if (k < r) {
            c = c * (r - k) / (k + 1);
        }
    }
    return cnt;
}

static long long powmod_i128(long long base, long long exp, long long mod) {
    i128 r = 1 % mod;
    i128 b = base % mod;
    while (exp) {
        if (exp & 1) r = r * b % mod;
        b = b * b % mod;
        exp >>= 1;
    }
    return (long long)r;
}

/* Compute A(n) given positions of 1-bits in n (sorted increasing), mod MOD */
static long long A_from_positions(long long *pos, int m, long long mod) {
    if (m == 0) return 0;
    if (m == 1) return 1 % mod;

    /* Precompute v_k: v_0=1, v_{k+1} = 5*v_k + 3, up to k = m-1 */
    long long *v = (long long*)malloc(m * sizeof(long long));
    v[0] = 1 % mod;
    for (int k = 1; k < m; k++) {
        v[k] = (5 * v[k-1] + 3) % mod;
    }

    long long ans = 1 % mod;
    /* desc = pos reversed (MSB -> LSB) */
    for (int i = 0; i < m - 1; i++) {
        /* gap between desc[i] and desc[i+1] */
        long long gap = pos[m - 1 - i] - pos[m - 1 - i - 1] - 1;
        if (gap <= 0) continue;
        long long base = v[i + 1];
        ans = (ans * powmod_i128(base, gap, mod)) % mod;
    }

    free(v);
    return ans;
}

long long p811_native(void) {
    long long t = 100000000000031LL; /* 10^14 + 31 */
    int r = 62;

    /* t is huge, so t >= max_binom_bitlen(r) is guaranteed.
       Total number of 1-bits: sum of popcount(C(62,k)) for k=0..62.
       Upper bound: 63 * 62 = 3906. Allocate generously. */
    int maxbits = 4000;
    long long *pos = (long long*)malloc(maxbits * sizeof(long long));
    int npos = one_positions_via_binom(t, r, pos);

    /* positions are already in increasing order from the binomial method */
    long long result = A_from_positions(pos, npos, MOD);

    free(pos);
    return result;
}
