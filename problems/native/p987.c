// Project Euler 987: Straights
// Count disjoint straights from one deck, excluding straight flushes.
// Answer for target=8 exceeds LLONG_MAX, so use unsigned + print function.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef __int128 i128;
typedef unsigned long long u64;

// 10 rank windows for straights. Ranks: 0=A,1=2,...,8=9,9=10,10=J,11=Q,12=K
static const int WINDOWS[10][5] = {
    {0,1,2,3,4}, {1,2,3,4,5}, {2,3,4,5,6}, {3,4,5,6,7}, {4,5,6,7,8},
    {5,6,7,8,9}, {6,7,8,9,10}, {7,8,9,10,11}, {8,9,10,11,12}, {9,10,11,12,0}
};

static int overlap[10][10];

static void build_overlap() {
    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++) {
            int shared = 0;
            for (int a = 0; a < 5 && !shared; a++)
                for (int b = 0; b < 5; b++)
                    if (WINDOWS[i][a] == WINDOWS[j][b]) { shared = 1; break; }
            overlap[i][j] = shared;
        }
}

// PERMS[n][k] = P(n,k) = n!/(n-k)!
static i128 PERMS[5][5];
static void build_perms() {
    for (int n = 0; n < 5; n++) {
        PERMS[n][0] = 1;
        i128 v = 1;
        for (int k = 1; k < 5; k++) {
            if (k <= n) { v *= (n - (k-1)); PERMS[n][k] = v; }
            else PERMS[n][k] = 0;
        }
    }
}

static i128 factorial(int n) {
    i128 r = 1;
    for (int i = 2; i <= n; i++) r *= i;
    return r;
}

// Compute colorings of all subsets for a set of labeled straights.
// dp[mask] = number of proper 4-colorings of the subgraph induced by mask.
static i128 coloring_dp[256];

static void colorings_of_all_subsets(const int *starts, int k) {
    int full = 1 << k;
    int adjacency[8] = {0};
    for (int i = 0; i < k; i++)
        for (int j = i+1; j < k; j++)
            if (overlap[starts[i]][starts[j]]) {
                adjacency[i] |= (1 << j);
                adjacency[j] |= (1 << i);
            }

    // independent[mask]
    char independent[256];
    independent[0] = 1;
    for (int mask = 1; mask < full; mask++) {
        int bit = mask & (-mask);
        int vertex = __builtin_ctz(bit);
        int rest = mask ^ bit;
        independent[mask] = independent[rest] && !(adjacency[vertex] & rest);
    }

    // DP over 4 colors
    i128 dp[256], new_dp[256];
    memset(dp, 0, sizeof(dp));
    dp[0] = 1;
    for (int c = 0; c < 4; c++) {
        for (int mask = 0; mask < full; mask++) {
            i128 total = 0;
            int sub = mask;
            while (1) {
                if (independent[sub])
                    total += dp[mask ^ sub];
                if (sub == 0) break;
                sub = (sub - 1) & mask;
            }
            new_dp[mask] = total;
        }
        memcpy(dp, new_dp, full * sizeof(i128));
    }
    memcpy(coloring_dp, dp, full * sizeof(i128));
}

// Count disjoint choices for labeled straight types, excluding straight flushes.
static i128 labeled_count(const int *starts, int k) {
    int full = 1 << k;
    colorings_of_all_subsets(starts, k);

    int total_active[13] = {0};
    int active_masks_by_rank[13] = {0};
    for (int idx = 0; idx < k; idx++) {
        int bit = 1 << idx;
        for (int r = 0; r < 5; r++) {
            int rank = WINDOWS[starts[idx]][r];
            total_active[rank]++;
            active_masks_by_rank[rank] |= bit;
        }
    }

    i128 total = 0;
    for (int mask = 0; mask < full; mask++) {
        i128 ways = 1;
        int ok = 1;
        for (int rank = 0; rank < 13 && ok; rank++) {
            int mono = __builtin_popcount(active_masks_by_rank[rank] & mask);
            int flex = total_active[rank] - mono;
            i128 wr = PERMS[4 - mono][flex];
            if (wr == 0) { ok = 0; break; }
            ways *= wr;
        }
        if (!ok) continue;
        i128 term = coloring_dp[mask] * ways;
        if (__builtin_popcount(mask) & 1)
            total -= term;
        else
            total += term;
    }
    return total;
}

// Backtracking for feasible type counts.
static int bt_counts[10];
static int bt_coverage[13];

// Callback for each feasible type count.
static i128 (*bt_callback)(const int *counts);
static i128 bt_result;

static void backtrack(int pos, int remaining) {
    if (pos == 10) {
        if (remaining == 0) {
            bt_result += bt_callback(bt_counts);
        }
        return;
    }
    for (int amount = 0; amount <= remaining; amount++) {
        int ok = 1;
        for (int r = 0; r < 5; r++) {
            bt_coverage[WINDOWS[pos][r]] += amount;
            if (bt_coverage[WINDOWS[pos][r]] > 4) ok = 0;
        }
        bt_counts[pos] = amount;
        if (ok)
            backtrack(pos + 1, remaining - amount);
        for (int r = 0; r < 5; r++)
            bt_coverage[WINDOWS[pos][r]] -= amount;
    }
    bt_counts[pos] = 0;
}

static i128 process_type_count(const int *counts) {
    int starts[8];
    int k = 0;
    i128 divisor = 1;
    for (int s = 0; s < 10; s++) {
        for (int a = 0; a < counts[s]; a++)
            starts[k++] = s;
        divisor *= factorial(counts[s]);
    }
    return labeled_count(starts, k) / divisor;
}

static i128 count_disjoint_straights(int target) {
    bt_callback = process_type_count;
    bt_result = 0;
    memset(bt_coverage, 0, sizeof(bt_coverage));
    backtrack(0, target);
    return bt_result;
}

void p987_print(void) {
    build_overlap();
    build_perms();
    i128 result = count_disjoint_straights(8);
    // Print as unsigned since it exceeds LLONG_MAX
    u64 val = (u64)result;
    printf("%llu\n", val);
}
