// Project Euler 879: Touch-screen Password
// Count distinct passwords on an n x n grid.
// Bitmask DP: dp[cur][mask] = number of non-empty continuations.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int gcd(int a, int b) {
    while (b) { int t = a % b; a = b; b = t; }
    return a;
}

static int between[16][16];

static void precompute_between(int n) {
    int N = n * n;
    memset(between, 0, sizeof(between));
    for (int a = 0; a < N; a++) {
        int xa = a % n, ya = a / n;
        for (int b = 0; b < N; b++) {
            if (a == b) continue;
            int xb = b % n, yb = b / n;
            int dx = xb - xa, dy = yb - ya;
            int g = gcd(dx < 0 ? -dx : dx, dy < 0 ? -dy : dy);
            if (g <= 1) continue;
            int sx = dx / g, sy = dy / g;
            int mask = 0;
            for (int k = 1; k < g; k++) {
                int x = xa + sx * k;
                int y = ya + sy * k;
                mask |= 1 << (y * n + x);
            }
            between[a][b] = mask;
        }
    }
}

// dp[cur][mask] - we use a flat array: dp[cur * 65536 + mask]
static long long dp[16 * 65536];

// buckets by popcount
static int buckets[17][65536];
static int bucket_count[17];

static int popcount(int x) {
    return __builtin_popcount(x);
}

static long long count_passwords(int n) {
    int N = n * n;
    if (N < 2) return 0;

    precompute_between(n);
    int all_mask = (1 << N) - 1;

    // Build buckets by popcount
    memset(bucket_count, 0, sizeof(bucket_count));
    for (int mask = 0; mask < (1 << N); mask++) {
        int pc = popcount(mask);
        buckets[pc][bucket_count[pc]++] = mask;
    }

    memset(dp, 0, sizeof(long long) * N * (1 << N));

    for (int k = N; k >= 1; k--) {
        for (int bi = 0; bi < bucket_count[k]; bi++) {
            int mask = buckets[k][bi];
            int remaining = all_mask ^ mask;
            if (remaining == 0) continue;

            int m = mask;
            while (m) {
                int lsb = m & (-m);
                int cur = __builtin_ctz(lsb);
                m ^= lsb;

                long long total = 0;
                int rem = remaining;
                while (rem) {
                    int bit = rem & (-rem);
                    int nxt = __builtin_ctz(bit);
                    rem ^= bit;

                    // nxt is directly selectable iff all intermediate points are already visited
                    if ((between[cur][nxt] & remaining) == 0) {
                        total += 1 + dp[nxt * (1 << N) + (mask | bit)];
                    }
                }
                dp[cur * (1 << N) + mask] = total;
            }
        }
    }

    long long total_passwords = 0;
    for (int start = 0; start < N; start++) {
        total_passwords += dp[start * (1 << N) + (1 << start)];
    }
    return total_passwords;
}

long long p879_native(void) {
    return count_passwords(4);
}
