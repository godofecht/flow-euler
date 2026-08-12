#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

/* Project Euler 917
   N x N matrix M[i,j] = a_i + b_j, path sum from (1,1) to (N,N) with Right/Down.
   A(N) = sum(a_i) + sum(b_j) + D(N), where D(N) is shortest path in grid graph.
   Only lower convex hull vertices matter. We compute hulls of (i, a_i) and (i, b_i),
   then DP on compressed grid.
*/

#define MOD 998388889LL
#define S1  102022661LL
#define N   10000000

static i128 cross(i128 x1, i128 y1, i128 x2, i128 y2, i128 x3, i128 y3) {
    /* cross product of (x2-x1, y2-y1) and (x3-x2, y3-y2) */
    return (x2 - x1) * (y3 - y2) - (y2 - y1) * (x3 - x2);
}

long long p917_native(void) {
    i64 s = S1;
    i128 sum_a = 0, sum_b = 0;

    /* Hull stacks: we store x and y for each hull vertex */
    /* For N=10^7, hulls are small (O(log) or O(sqrt) typically) */
    i64 *ax = malloc(sizeof(i64) * 1024);
    i64 *ay = malloc(sizeof(i64) * 1024);
    i64 *bx = malloc(sizeof(i64) * 1024);
    i64 *by = malloc(sizeof(i64) * 1024);
    int asz = 0, bsz = 0;
    int acap = 1024, bcap = 1024;

    for (i64 i = 1; i <= N; i++) {
        i64 a = s;
        s = (s * s) % MOD;
        i64 b = s;
        s = (s * s) % MOD;

        sum_a += a;
        sum_b += b;

        /* Update lower hull for (i, a) */
        while (asz >= 2) {
            i128 cr = cross(ax[asz-2], ay[asz-2], ax[asz-1], ay[asz-1], i, a);
            if (cr <= 0) {
                asz--;
            } else {
                break;
            }
        }
        if (asz >= acap) {
            acap *= 2;
            ax = realloc(ax, sizeof(i64) * acap);
            ay = realloc(ay, sizeof(i64) * acap);
        }
        ax[asz] = i;
        ay[asz] = a;
        asz++;

        /* Update lower hull for (i, b) */
        while (bsz >= 2) {
            i128 cr = cross(bx[bsz-2], by[bsz-2], bx[bsz-1], by[bsz-1], i, b);
            if (cr <= 0) {
                bsz--;
            } else {
                break;
            }
        }
        if (bsz >= bcap) {
            bcap *= 2;
            bx = realloc(bx, sizeof(i64) * bcap);
            by = realloc(by, sizeof(i64) * bcap);
        }
        bx[bsz] = i;
        by[bsz] = b;
        bsz++;
    }

    /* DP on compressed grid */
    int R = asz, C = bsz;
    i128 *dp = malloc(sizeof(i128) * C);
    i128 *ndp = malloc(sizeof(i128) * C);

    dp[0] = 0;
    /* First row: only moves right, cost = a0 * (col distance) */
    i64 a0 = ay[0];
    for (int j = 1; j < C; j++) {
        dp[j] = dp[j-1] + (i128)a0 * (bx[j] - bx[j-1]);
    }

    for (int i = 1; i < R; i++) {
        i128 dr = (i128)(ax[i] - ax[i-1]);
        /* First column: only moves down */
        ndp[0] = dp[0] + (i128)by[0] * dr;

        i64 ai = ay[i];
        for (int j = 1; j < C; j++) {
            i128 down = dp[j] + (i128)by[j] * dr;
            i128 dc = (i128)(bx[j] - bx[j-1]);
            i128 right = ndp[j-1] + (i128)ai * dc;
            ndp[j] = down < right ? down : right;
        }
        /* swap dp and ndp */
        i128 *tmp = dp;
        dp = ndp;
        ndp = tmp;
    }

    i128 D = dp[C-1];
    i128 result = sum_a + sum_b + D;

    free(dp); free(ndp);
    free(ax); free(ay); free(bx); free(by);

    return (i64)result;
}
