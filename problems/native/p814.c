// Project Euler 814: Mezzo-forte
// Count assignments of looking directions for 4n people on a circle
// such that exactly half scream (exactly n mutual pairs).
// DP with width-2 transfer over 2n columns in a 2x(2n) twisted ladder.
// Port of the Python reference solver to C.
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define MOD 998244353LL

static int swap2bits(int x) {
    return ((x & 1) << 1) | ((x >> 1) & 1);
}

// Transition groups: for each incoming mask (0..3), list of (out_mask, delta, count)
typedef struct { int out_mask, delta, cnt; } Group;

static const Group GROUPS[4][7] = {
    // mask 0: in_top=0, in_bottom=0
    {{0,0,3},{2,0,2},{1,0,2},{3,0,1},{0,1,1},{-1,-1,-1},{-1,-1,-1}},
    // mask 1: in_top=1, in_bottom=0
    {{0,1,3},{2,1,1},{1,0,2},{3,0,1},{0,0,1},{2,0,1},{-1,-1,-1}},
    // mask 2: in_top=0, in_bottom=1
    {{0,1,3},{2,0,2},{0,0,1},{1,1,1},{3,0,1},{1,0,1},{-1,-1,-1}},
    // mask 3: in_top=1, in_bottom=1
    {{0,2,1},{2,1,1},{0,1,3},{1,1,1},{3,0,1},{1,0,1},{2,0,1}},
};

static int group_sizes[4] = {5, 6, 6, 7};

long long p814_native(void) {
    int n = 1000;
    int target = n;
    int cols = 2 * n;

    // Build GROUPS_LAST (swap out_mask bits)
    Group GROUPS_LAST[4][7];
    for (int mi = 0; mi < 4; mi++) {
        for (int j = 0; j < group_sizes[mi]; j++) {
            GROUPS_LAST[mi][j].out_mask = swap2bits(GROUPS[mi][j].out_mask);
            GROUPS_LAST[mi][j].delta = GROUPS[mi][j].delta;
            GROUPS_LAST[mi][j].cnt = GROUPS[mi][j].cnt;
        }
    }

    long long total = 0;

    // dp[mask][k] and dp_next[mask][k]
    long long *dp = malloc(4 * (target + 1) * sizeof(long long));
    long long *dpn = malloc(4 * (target + 1) * sizeof(long long));
    if (!dp || !dpn) return 0;

    for (int init_in = 0; init_in < 4; init_in++) {
        memset(dp, 0, 4 * (target + 1) * sizeof(long long));
        dp[init_in * (target + 1) + 0] = 1;

        // Normal columns 0 .. cols-2
        for (int col = 0; col < cols - 1; col++) {
            memset(dpn, 0, 4 * (target + 1) * sizeof(long long));
            for (int mask_in = 0; mask_in < 4; mask_in++) {
                long long *arr = &dp[mask_in * (target + 1)];
                for (int gi = 0; gi < group_sizes[mask_in]; gi++) {
                    int out_mask = GROUPS[mask_in][gi].out_mask;
                    int delta = GROUPS[mask_in][gi].delta;
                    int cnt = GROUPS[mask_in][gi].cnt;
                    long long *dest = &dpn[out_mask * (target + 1)];

                    if (cnt == 1) {
                        if (delta == 0) {
                            for (int k = 0; k <= target; k++) {
                                long long x = dest[k] + arr[k];
                                if (x >= MOD) x -= MOD;
                                dest[k] = x;
                            }
                        } else if (delta == 1) {
                            for (int k = 0; k < target; k++) {
                                int idx = k + 1;
                                long long x = dest[idx] + arr[k];
                                if (x >= MOD) x -= MOD;
                                dest[idx] = x;
                            }
                        } else {
                            for (int k = 0; k < target - 1; k++) {
                                int idx = k + 2;
                                long long x = dest[idx] + arr[k];
                                if (x >= MOD) x -= MOD;
                                dest[idx] = x;
                            }
                        }
                    } else if (cnt == 2) {
                        if (delta == 0) {
                            for (int k = 0; k <= target; k++) {
                                long long add = arr[k] * 2;
                                if (add >= MOD) add -= MOD;
                                long long x = dest[k] + add;
                                if (x >= MOD) x -= MOD;
                                dest[k] = x;
                            }
                        } else if (delta == 1) {
                            for (int k = 0; k < target; k++) {
                                long long add = arr[k] * 2;
                                if (add >= MOD) add -= MOD;
                                int idx = k + 1;
                                long long x = dest[idx] + add;
                                if (x >= MOD) x -= MOD;
                                dest[idx] = x;
                            }
                        } else {
                            for (int k = 0; k < target - 1; k++) {
                                long long add = arr[k] * 2;
                                if (add >= MOD) add -= MOD;
                                int idx = k + 2;
                                long long x = dest[idx] + add;
                                if (x >= MOD) x -= MOD;
                                dest[idx] = x;
                            }
                        }
                    } else { // cnt == 3
                        if (delta == 0) {
                            for (int k = 0; k <= target; k++) {
                                long long add = arr[k] * 3;
                                if (add >= MOD) add -= MOD;
                                if (add >= MOD) add -= MOD;
                                long long x = dest[k] + add;
                                if (x >= MOD) x -= MOD;
                                dest[k] = x;
                            }
                        } else if (delta == 1) {
                            for (int k = 0; k < target; k++) {
                                long long add = arr[k] * 3;
                                if (add >= MOD) add -= MOD;
                                if (add >= MOD) add -= MOD;
                                int idx = k + 1;
                                long long x = dest[idx] + add;
                                if (x >= MOD) x -= MOD;
                                dest[idx] = x;
                            }
                        } else {
                            for (int k = 0; k < target - 1; k++) {
                                long long add = arr[k] * 3;
                                if (add >= MOD) add -= MOD;
                                if (add >= MOD) add -= MOD;
                                int idx = k + 2;
                                long long x = dest[idx] + add;
                                if (x >= MOD) x -= MOD;
                                dest[idx] = x;
                            }
                        }
                    }
                }
            }
            long long *tmp = dp; dp = dpn; dpn = tmp;
        }

        // Final column (wrap-around with row swap)
        memset(dpn, 0, 4 * (target + 1) * sizeof(long long));
        for (int mask_in = 0; mask_in < 4; mask_in++) {
            long long *arr = &dp[mask_in * (target + 1)];
            for (int gi = 0; gi < group_sizes[mask_in]; gi++) {
                int next_in = GROUPS_LAST[mask_in][gi].out_mask;
                int delta = GROUPS_LAST[mask_in][gi].delta;
                int cnt = GROUPS_LAST[mask_in][gi].cnt;
                long long *dest = &dpn[next_in * (target + 1)];

                if (cnt == 1) {
                    if (delta == 0) {
                        for (int k = 0; k <= target; k++) {
                            long long x = dest[k] + arr[k];
                            if (x >= MOD) x -= MOD;
                            dest[k] = x;
                        }
                    } else if (delta == 1) {
                        for (int k = 0; k < target; k++) {
                            int idx = k + 1;
                            long long x = dest[idx] + arr[k];
                            if (x >= MOD) x -= MOD;
                            dest[idx] = x;
                        }
                    } else {
                        for (int k = 0; k < target - 1; k++) {
                            int idx = k + 2;
                            long long x = dest[idx] + arr[k];
                            if (x >= MOD) x -= MOD;
                            dest[idx] = x;
                        }
                    }
                } else if (cnt == 2) {
                    if (delta == 0) {
                        for (int k = 0; k <= target; k++) {
                            long long add = arr[k] * 2;
                            if (add >= MOD) add -= MOD;
                            long long x = dest[k] + add;
                            if (x >= MOD) x -= MOD;
                            dest[k] = x;
                        }
                    } else if (delta == 1) {
                        for (int k = 0; k < target; k++) {
                            long long add = arr[k] * 2;
                            if (add >= MOD) add -= MOD;
                            int idx = k + 1;
                            long long x = dest[idx] + add;
                            if (x >= MOD) x -= MOD;
                            dest[idx] = x;
                        }
                    } else {
                        for (int k = 0; k < target - 1; k++) {
                            long long add = arr[k] * 2;
                            if (add >= MOD) add -= MOD;
                            int idx = k + 2;
                            long long x = dest[idx] + add;
                            if (x >= MOD) x -= MOD;
                            dest[idx] = x;
                        }
                    }
                } else { // cnt == 3
                    if (delta == 0) {
                        for (int k = 0; k <= target; k++) {
                            long long add = arr[k] * 3;
                            if (add >= MOD) add -= MOD;
                            if (add >= MOD) add -= MOD;
                            long long x = dest[k] + add;
                            if (x >= MOD) x -= MOD;
                            dest[k] = x;
                        }
                    } else if (delta == 1) {
                        for (int k = 0; k < target; k++) {
                            long long add = arr[k] * 3;
                            if (add >= MOD) add -= MOD;
                            if (add >= MOD) add -= MOD;
                            int idx = k + 1;
                            long long x = dest[idx] + add;
                            if (x >= MOD) x -= MOD;
                            dest[idx] = x;
                        }
                    } else {
                        for (int k = 0; k < target - 1; k++) {
                            long long add = arr[k] * 3;
                            if (add >= MOD) add -= MOD;
                            if (add >= MOD) add -= MOD;
                            int idx = k + 2;
                            long long x = dest[idx] + add;
                            if (x >= MOD) x -= MOD;
                            dest[idx] = x;
                        }
                    }
                }
            }
        }
        long long *tmp = dp; dp = dpn; dpn = tmp;

        total += dp[init_in * (target + 1) + target];
        total %= MOD;
    }

    free(dp);
    free(dpn);
    return total;
}
