/* Project Euler 859: Cookie Game
 *
 * Odd (Left) and Even (Right) play on piles of cookies.
 * A move chooses a single pile and replaces it by two equal smaller piles
 * (after eating 1 cookie for Odd, or 2 cookies for Even).
 *
 * We count C(N): the number of initial (unordered) pile partitions of N
 * for which Even has a winning strategy (Odd starts).
 *
 * Each pile is a cold combinatorial game whose value g(n) is an integer.
 * The total game value is the sum of pile values.  Even (Right) wins iff
 * the total value is <= 0 (value 0 means the second player wins).
 *
 * A coin-change style DP counts unordered partitions of N grouped by
 * total value, then sums the counts with value <= 0.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

long long p859_native(void) {
    int n = 300;

    /* g[k]: combinatorial game value of a single pile of size k. */
    long long *g = calloc((size_t)(n + 1), sizeof(long long));
    g[0] = 0;
    for (int k = 1; k <= n; k++) {
        if (k & 1) {
            /* k = 2m + 1: Left moves to two piles of size m. */
            int m = k >> 1;
            long long option = 2 * g[m];
            g[k] = (option < 0) ? 0 : option + 1;
        } else {
            /* k = 2m: Right moves to two piles of size m - 1. */
            int m = k >> 1;
            long long option = 2 * g[m - 1];
            g[k] = (option > 0) ? 0 : option - 1;
        }
    }

    /* Value range for the full position: [-n/2, n]. */
    int v_min = -(n / 2);
    int v_max = n;
    int width = v_max - v_min + 1;
    int offset = -v_min; /* index where value 0 lives */

    /* dp[cookies * width + vi] = #partitions of `cookies` with total value
     * mapped to index vi. */
    long long *dp = calloc((size_t)(n + 1) * (size_t)width, sizeof(long long));
    dp[offset] = 1; /* dp[0]: empty partition, value 0. */

    for (int size = 1; size <= n; size++) {
        long long dv = g[size];
        if (dv >= 0) {
            int d = (int)dv;
            for (int cookies = size; cookies <= n; cookies++) {
                long long *src = dp + (long long)(cookies - size) * width;
                long long *dst = dp + (long long)cookies * width;
                for (int vi = d; vi < width; vi++) {
                    dst[vi] += src[vi - d];
                }
            }
        } else {
            int shift = (int)(-dv);
            int limit = width - shift;
            for (int cookies = size; cookies <= n; cookies++) {
                long long *src = dp + (long long)(cookies - size) * width;
                long long *dst = dp + (long long)cookies * width;
                for (int vi = 0; vi < limit; vi++) {
                    dst[vi] += src[vi + shift];
                }
            }
        }
    }

    /* Sum counts with total value <= 0 (indices 0..offset). */
    long long total = 0;
    long long *row = dp + (long long)n * width;
    for (int vi = 0; vi <= offset; vi++) {
        total += row[vi];
    }

    free(g);
    free(dp);
    return total;
}
