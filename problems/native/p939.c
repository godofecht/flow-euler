/* Project Euler 939: Partisan Nim
 * Compute E(N) mod 1234567891 for N=5000.
 * O(N^2) partition DP.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef uint32_t u32;
typedef uint64_t u64;
typedef int64_t i64;

#define MOD 1234567891ULL
#define N 5000

static u32 *cnt0, *cnt1;
static u32 *start;

static void build_partition_stats(void) {
    /* Triangular offsets: start[n] = n*(n+1)/2 */
    start = (u32*)malloc((N + 1) * sizeof(u32));
    start[0] = 0;
    for (u32 n = 1; n <= N; n++) {
        start[n] = start[n - 1] + n;
    }

    u32 total_entries = (u32)((N + 1) * (N + 2) / 2);
    cnt0 = (u32*)calloc(total_entries, sizeof(u32));
    cnt1 = (u32*)calloc(total_entries, sizeof(u32));

    /* Empty partition of 0 */
    cnt0[0] = 1;

    u32 *prev0 = (u32*)calloc(N + 1, sizeof(u32));
    u32 *prev1 = (u32*)calloc(N + 1, sizeof(u32));
    prev0[0] = 1; /* f[0,0,0] = 1 */

    u32 *cur0 = (u32*)malloc((N + 1) * sizeof(u32));
    u32 *cur1 = (u32*)malloc((N + 1) * sizeof(u32));

    for (u32 k = 1; k <= N; k++) {
        memset(cur0, 0, (N + 1) * sizeof(u32));
        memset(cur1, 0, (N + 1) * sizeof(u32));

        if ((k & 1) == 0) {
            /* parity unchanged in (n-k, k) branch */
            for (u32 n = k; n <= N; n++) {
                u32 a0 = prev1[n - 1];
                u32 a1 = prev0[n - 1];
                u32 b0 = cur0[n - k];
                u32 b1 = cur1[n - k];

                u32 s0 = a0 + b0;
                if (s0 >= MOD) s0 -= MOD;
                u32 s1 = a1 + b1;
                if (s1 >= MOD) s1 -= MOD;

                cur0[n] = s0;
                cur1[n] = s1;

                u32 v = n - k;
                u32 idx = start[n] + v;

                u32 t = cnt0[idx] + s0;
                if (t >= MOD) t -= MOD;
                cnt0[idx] = t;

                t = cnt1[idx] + s1;
                if (t >= MOD) t -= MOD;
                cnt1[idx] = t;
            }
        } else {
            /* parity flips in (n-k, k) branch */
            for (u32 n = k; n <= N; n++) {
                u32 a0 = prev1[n - 1];
                u32 a1 = prev0[n - 1];
                u32 b0 = cur1[n - k];
                u32 b1 = cur0[n - k];

                u32 s0 = a0 + b0;
                if (s0 >= MOD) s0 -= MOD;
                u32 s1 = a1 + b1;
                if (s1 >= MOD) s1 -= MOD;

                cur0[n] = s0;
                cur1[n] = s1;

                u32 v = n - k;
                u32 idx = start[n] + v;

                u32 t = cnt0[idx] + s0;
                if (t >= MOD) t -= MOD;
                cnt0[idx] = t;

                t = cnt1[idx] + s1;
                if (t >= MOD) t -= MOD;
                cnt1[idx] = t;
            }
        }

        /* swap prev and cur */
        u32 *tmp = prev0; prev0 = cur0; cur0 = tmp;
        tmp = prev1; prev1 = cur1; cur1 = tmp;
    }

    free(prev0); free(prev1);
    free(cur0); free(cur1);
}

static u64 compute_E(void) {
    build_partition_stats();

    u32 *cum0 = (u32*)calloc(N + 1, sizeof(u32));
    u32 *cum1 = (u32*)calloc(N + 1, sizeof(u32));
    u32 *pref0 = (u32*)malloc((N + 1) * sizeof(u32));
    u32 *pref1 = (u32*)malloc((N + 1) * sizeof(u32));

    u64 ans = 0;

    for (u32 m = 0; m <= N; m++) {
        u32 sb = start[m];

        /* Add distributions for total b = m into cumulative B-side counts */
        for (u32 v = 0; v <= m; v++) {
            u32 idx = sb + v;

            u32 x = cnt0[idx];
            if (x) {
                u32 s = cum0[v] + x;
                if (s >= MOD) s -= MOD;
                cum0[v] = s;
            }

            x = cnt1[idx];
            if (x) {
                u32 s = cum1[v] + x;
                if (s >= MOD) s -= MOD;
                cum1[v] = s;
            }
        }

        u32 a = N - m;
        u32 sa = start[a];

        /* Build prefix sums up to v = a */
        u32 run = 0;
        for (u32 v = 0; v <= a; v++) {
            run += cum0[v];
            if (run >= MOD) run -= MOD;
            pref0[v] = run;
        }

        run = 0;
        for (u32 v = 0; v <= a; v++) {
            run += cum1[v];
            if (run >= MOD) run -= MOD;
            pref1[v] = run;
        }

        /* Count A partitions of total a against all B partitions of total <= m */
        for (u32 vA = 0; vA <= a; vA++) {
            u32 idxA = sa + vA;
            u32 a0 = cnt0[idxA];
            u32 a1 = cnt1[idxA];

            int x_same = (int)vA - 1;
            int x_diff = (int)vA - 2;

            if (a0) {
                u64 s = 0;
                if (x_same >= 0) s += pref0[x_same];
                if (x_diff >= 0) s += pref1[x_diff];
                s %= MOD;
                ans = (ans + (u64)a0 * s) % MOD;
            }

            if (a1) {
                u64 t = 0;
                if (x_same >= 0) t += pref1[x_same];
                if (x_diff >= 0) t += pref0[x_diff];
                t %= MOD;
                ans = (ans + (u64)a1 * t) % MOD;
            }
        }
    }

    free(cum0); free(cum1);
    free(pref0); free(pref1);
    free(cnt0); free(cnt1);
    free(start);

    return ans;
}

long long p939_native(void) {
    return (long long)compute_E();
}
