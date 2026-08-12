/* Project Euler 928: Cribbage
 *
 * Count non-empty hands (subsets of 52-card deck) where
 * hand_score == cribbage_score.
 *
 * Meet-in-the-middle: left ranks 1..7, right ranks 8..13.
 * Left: 5^7 = 78125 combos. Right: 5^6 = 15625 combos.
 *
 * Condition: A_L + A_R - run_boundary == 2 * dp15_total
 * where dp15_total = L15 + L7*c8 + L6*c9 + L5*n10
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;

static const int COMB4[5] = {1, 4, 6, 4, 1};
static int CHOOSE[5][5];

static void init_choose(void) {
    for (int n = 0; n < 5; n++) {
        CHOOSE[n][0] = 1;
        for (int k = 1; k <= n; k++)
            CHOOSE[n][k] = CHOOSE[n-1][k-1] + (k <= n-1 ? CHOOSE[n-1][k] : 0);
    }
}

static int rank_value(int r) { return r >= 10 ? 10 : r; }

/* Subset-sum DP update: multiply by sum_{k=0..c} C(c,k) x^(k*v) */
static void dp_update(i64 *dp, int v, int c) {
    if (c == 0) return;
    i64 nd[16] = {0};
    for (int s = 0; s < 16; s++) {
        i64 tot = 0;
        int maxk = c;
        if (v > 0 && s / v < maxk) maxk = s / v;
        for (int k = 0; k <= maxk; k++)
            tot += dp[s - k * v] * CHOOSE[c][k];
        nd[s] = tot;
    }
    memcpy(dp, nd, 16 * sizeof(i64));
}

/* ---- Left entry ---- */
typedef struct {
    int L5, L6, L7, L15;
    int tail_len, tail_prod;
    int A_L;
    i64 weight;
} LeftEntry;

/* ---- Right entry ---- */
typedef struct {
    int c8, c9, n10;
    int head_len, head_prod;
    int A_R;
    i64 weight;
} RightEntry;

static LeftEntry left_entries[78125];
static int n_left = 0;
static RightEntry right_entries[15625];
static int n_right = 0;

static int cmp_left(const void *a, const void *b) {
    const LeftEntry *x = a, *y = b;
    if (x->L5 != y->L5) return x->L5 - y->L5;
    if (x->L6 != y->L6) return x->L6 - y->L6;
    if (x->L7 != y->L7) return x->L7 - y->L7;
    if (x->L15 != y->L15) return x->L15 - y->L15;
    if (x->tail_len != y->tail_len) return x->tail_len - y->tail_len;
    if (x->tail_prod != y->tail_prod) return x->tail_prod - y->tail_prod;
    return x->A_L - y->A_L;
}

static int cmp_right(const void *a, const void *b) {
    const RightEntry *x = a, *y = b;
    if (x->c8 != y->c8) return x->c8 - y->c8;
    if (x->c9 != y->c9) return x->c9 - y->c9;
    if (x->n10 != y->n10) return x->n10 - y->n10;
    if (x->head_len != y->head_len) return x->head_len - y->head_len;
    if (x->head_prod != y->head_prod) return x->head_prod - y->head_prod;
    return x->A_R - y->A_R;
}

static int run_boundary(int tl, int tp, int hl, int hp) {
    int L, P;
    if (tl > 0 && hl > 0) { L = tl + hl; P = tp * hp; }
    else if (tl > 0) { L = tl; P = tp; }
    else if (hl > 0) { L = hl; P = hp; }
    else return 0;
    return (L >= 3) ? L * P : 0;
}

static void enumerate_left(void) {
    for (int c1 = 0; c1 < 5; c1++)
    for (int c2 = 0; c2 < 5; c2++)
    for (int c3 = 0; c3 < 5; c3++)
    for (int c4 = 0; c4 < 5; c4++)
    for (int c5 = 0; c5 < 5; c5++)
    for (int c6 = 0; c6 < 5; c6++)
    for (int c7 = 0; c7 < 5; c7++) {
        i64 w = (i64)COMB4[c1] * COMB4[c2] * COMB4[c3] * COMB4[c4]
                * COMB4[c5] * COMB4[c6] * COMB4[c7];
        if (w == 0) continue;

        int counts[7] = {c1, c2, c3, c4, c5, c6, c7};

        int hv = 1*c1 + 2*c2 + 3*c3 + 4*c4 + 5*c5 + 6*c6 + 7*c7;
        int pairs = c1*(c1-1) + c2*(c2-1) + c3*(c3-1) + c4*(c4-1)
                  + c5*(c5-1) + c6*(c6-1) + c7*(c7-1);

        i64 dp[16] = {0};
        dp[0] = 1;
        dp_update(dp, 1, c1);
        dp_update(dp, 2, c2);
        dp_update(dp, 3, c3);
        dp_update(dp, 4, c4);
        dp_update(dp, 5, c5);
        dp_update(dp, 6, c6);
        dp_update(dp, 7, c7);

        int L5 = (int)dp[5], L6 = (int)dp[6], L7 = (int)dp[7], L15 = (int)dp[15];
        if (L15 > 170) continue;

        /* run_internal: exclude trailing segment reaching rank 7 */
        int run_internal = 0, cur_len = 0, cur_prod = 1;
        for (int i = 0; i < 7; i++) {
            if (counts[i] == 0) {
                if (cur_len >= 3) run_internal += cur_len * cur_prod;
                cur_len = 0; cur_prod = 1;
            } else {
                cur_len++;
                cur_prod *= counts[i];
            }
        }
        int tail_len = cur_len;
        int tail_prod = cur_len > 0 ? cur_prod : 1;

        int A_L = hv - pairs - run_internal;

        left_entries[n_left].L5 = L5;
        left_entries[n_left].L6 = L6;
        left_entries[n_left].L7 = L7;
        left_entries[n_left].L15 = L15;
        left_entries[n_left].tail_len = tail_len;
        left_entries[n_left].tail_prod = tail_prod;
        left_entries[n_left].A_L = A_L;
        left_entries[n_left].weight = w;
        n_left++;
    }
}

static void enumerate_right(void) {
    for (int c8 = 0; c8 < 5; c8++)
    for (int c9 = 0; c9 < 5; c9++)
    for (int c10 = 0; c10 < 5; c10++)
    for (int c11 = 0; c11 < 5; c11++)
    for (int c12 = 0; c12 < 5; c12++)
    for (int c13 = 0; c13 < 5; c13++) {
        i64 w = (i64)COMB4[c8] * COMB4[c9] * COMB4[c10] * COMB4[c11]
                * COMB4[c12] * COMB4[c13];
        if (w == 0) continue;

        int n10 = c10 + c11 + c12 + c13;
        int hv = 8*c8 + 9*c9 + 10*n10;
        int pairs = c8*(c8-1) + c9*(c9-1) + c10*(c10-1) + c11*(c11-1)
                  + c12*(c12-1) + c13*(c13-1);

        int counts[6] = {c8, c9, c10, c11, c12, c13};

        /* head segment starting at rank 8 */
        int head_len = 0, head_prod = 1;
        for (int i = 0; i < 6; i++) {
            if (counts[i] == 0) break;
            head_len++;
            head_prod *= counts[i];
        }
        if (head_len == 0) { head_len = 0; head_prod = 1; }

        /* run_internal: all segments length>=3 EXCEPT the one starting at rank 8 */
        int run_internal = 0, cur_len = 0, cur_prod = 1, seg_start = -1;
        for (int i = 0; i < 6; i++) {
            if (counts[i] == 0) {
                if (cur_len >= 3 && seg_start != 0)
                    run_internal += cur_len * cur_prod;
                cur_len = 0; cur_prod = 1; seg_start = -1;
            } else {
                if (cur_len == 0) seg_start = i;
                cur_len++;
                cur_prod *= counts[i];
            }
        }
        if (cur_len >= 3 && seg_start != 0)
            run_internal += cur_len * cur_prod;

        int A_R = hv - pairs - run_internal;

        right_entries[n_right].c8 = c8;
        right_entries[n_right].c9 = c9;
        right_entries[n_right].n10 = n10;
        right_entries[n_right].head_len = head_len;
        right_entries[n_right].head_prod = head_prod;
        right_entries[n_right].A_R = A_R;
        right_entries[n_right].weight = w;
        n_right++;
    }
}

/* Binary search for A_R value in a sorted range of right entries */
static i64 binary_search_right(int lo, int hi, int target_A_R) {
    /* right_entries[lo..hi-1] are sorted by A_R (within same c8,c9,n10,head_len,head_prod) */
    int orig_hi = hi;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (right_entries[mid].A_R < target_A_R) lo = mid + 1;
        else hi = mid;
    }
    if (lo < orig_hi && right_entries[lo].A_R == target_A_R)
        return right_entries[lo].weight;
    return 0;
}

long long p928_native(void) {
    init_choose();
    enumerate_left();
    enumerate_right();

    /* Sort both arrays */
    qsort(left_entries, n_left, sizeof(LeftEntry), cmp_left);
    qsort(right_entries, n_right, sizeof(RightEntry), cmp_right);

    /* Merge left entries with same key and A_L */
    int left_group_start[78125];
    int left_group_end[78125];
    int n_left_groups = 0;

    int i = 0;
    while (i < n_left) {
        int j = i;
        while (j < n_left && cmp_left(&left_entries[i], &left_entries[j]) == 0) j++;
        /* Merge weights for same (key, A_L) */
        i64 total_w = 0;
        for (int k = i; k < j; k++) total_w += left_entries[k].weight;
        left_entries[i].weight = total_w;
        left_group_start[n_left_groups] = i;
        left_group_end[n_left_groups] = j;
        n_left_groups++;
        i = j;
    }

    /* Build right group index: for each unique (c8, c9, n10, head_len, head_prod),
       store start and end in the sorted array (entries within are sorted by A_R) */
    typedef struct { int c8, c9, n10, head_len, head_prod, start, end; } RightGroup;
    RightGroup right_groups[15625];
    int n_right_groups = 0;

    i = 0;
    while (i < n_right) {
        int j = i;
        while (j < n_right &&
               right_entries[j].c8 == right_entries[i].c8 &&
               right_entries[j].c9 == right_entries[i].c9 &&
               right_entries[j].n10 == right_entries[i].n10 &&
               right_entries[j].head_len == right_entries[i].head_len &&
               right_entries[j].head_prod == right_entries[i].head_prod)
            j++;
        /* Merge weights for same A_R within this group */
        int k = i;
        while (k < j) {
            int l = k;
            while (l < j && right_entries[l].A_R == right_entries[k].A_R) l++;
            i64 total_w = 0;
            for (int m = k; m < l; m++) total_w += right_entries[m].weight;
            right_entries[k].weight = total_w;
            k = l;
        }
        right_groups[n_right_groups].c8 = right_entries[i].c8;
        right_groups[n_right_groups].c9 = right_entries[i].c9;
        right_groups[n_right_groups].n10 = right_entries[i].n10;
        right_groups[n_right_groups].head_len = right_entries[i].head_len;
        right_groups[n_right_groups].head_prod = right_entries[i].head_prod;
        right_groups[n_right_groups].start = i;
        right_groups[n_right_groups].end = j;
        n_right_groups++;
        i = j;
    }

    /* Build right (c8, c9, n10) super-group index */
    typedef struct { int c8, c9, n10, start, end; } RightSuperGroup;
    RightSuperGroup right_super[500];
    int n_right_super = 0;

    i = 0;
    while (i < n_right_groups) {
        int j = i;
        while (j < n_right_groups &&
               right_groups[j].c8 == right_groups[i].c8 &&
               right_groups[j].c9 == right_groups[i].c9 &&
               right_groups[j].n10 == right_groups[i].n10)
            j++;
        right_super[n_right_super].c8 = right_groups[i].c8;
        right_super[n_right_super].c9 = right_groups[i].c9;
        right_super[n_right_super].n10 = right_groups[i].n10;
        right_super[n_right_super].start = i;
        right_super[n_right_super].end = j;
        n_right_super++;
        i = j;
    }

    /* Combination step */
    i64 total = 0;

    for (int lg = 0; lg < n_left_groups; lg++) {
        LeftEntry *le = &left_entries[left_group_start[lg]];
        int L5 = le->L5, L6 = le->L6, L7 = le->L7, L15 = le->L15;
        int tl = le->tail_len, tp = le->tail_prod;

        /* Left group has entries from left_group_start[lg] to left_group_end[lg]-1
           All have same A_L (since we sorted by key then A_L, and merged) */
        int A_L = le->A_L;
        i64 wL = le->weight;

        for (int rs = 0; rs < n_right_super; rs++) {
            int c8 = right_super[rs].c8;
            int c9 = right_super[rs].c9;
            int n10 = right_super[rs].n10;
            int dp15 = L15 + L7 * c8 + L6 * c9 + L5 * n10;
            if (dp15 > 170) continue;
            int base = 2 * dp15;

            for (int rg = right_super[rs].start; rg < right_super[rs].end; rg++) {
                int hl = right_groups[rg].head_len;
                int hp = right_groups[rg].head_prod;
                int req = base + run_boundary(tl, tp, hl, hp);
                int want = req - A_L;
                i64 wR = binary_search_right(right_groups[rg].start,
                                             right_groups[rg].end, want);
                if (wR) total += wL * wR;
            }
        }
    }

    /* Exclude empty hand */
    total -= 1;
    return total;
}
