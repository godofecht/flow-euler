#include <stdint.h>
#include <string.h>
#include <stdlib.h>

typedef long long i64;
typedef unsigned long long u64;

#define L 36

/* ---------- GCD, extended GCD ---------- */

static i64 gcd64(i64 a, i64 b) {
    while (b) { i64 t = b; b = a % b; a = t; }
    return a;
}

static void egcd(i64 a, i64 b, i64 *g, i64 *x, i64 *y) {
    i64 x0 = 1, y0 = 0, x1 = 0, y1 = 1;
    while (b != 0) {
        i64 q = a / b;
        i64 t = b; b = a % b; a = t;
        i64 tx = x1; x1 = x0 - q * x1; x0 = tx;
        i64 ty = y1; y1 = y0 - q * y1; y0 = ty;
    }
    *g = a; *x = x0; *y = y0;
}

/* ---------- CRT merge ---------- */

static int crt_merge(i64 r1, i64 m1, i64 r2, i64 m2, i64 *r_out, i64 *m_out) {
    i64 g = gcd64(m1, m2);
    i64 diff_raw = r2 - r1;
    if (diff_raw % g != 0)
        return 0;
    i64 l = (m1 / g) * m2;
    i64 m1g = m1 / g;
    i64 m2g = m2 / g;
    i64 diff = diff_raw / g;
    i64 eg, einv, ey;
    egcd(m1g, m2g, &eg, &einv, &ey);
    einv = ((einv % m2g) + m2g) % m2g;
    i64 t = ((diff * einv) % m2g + m2g) % m2g;
    i64 r = (r1 + m1 * t) % l;
    *r_out = r;
    *m_out = l;
    return 1;
}

/* ---------- lcm ---------- */

static i64 lcm_upto(int n) {
    i64 m = 1;
    for (int i = 1; i <= n; i++)
        m = (m / gcd64(m, i)) * i;
    return m;
}

/* ---------- Hash table for visited states ---------- */

#define HT_SIZE (1 << 18)
#define HT_MASK (HT_SIZE - 1)

typedef struct { i64 r, m, unused, remaining; } state_t;
static state_t g_ht[HT_SIZE];
static int g_ht_used[HT_SIZE];

static u64 state_hash(i64 r, i64 m, i64 unused, i64 remaining) {
    u64 h = (u64)r;
    h ^= (u64)m + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    h ^= (u64)unused + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    h ^= (u64)remaining + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
}

static int ht_contains(i64 r, i64 m, i64 unused, i64 remaining) {
    u64 h = state_hash(r, m, unused, remaining) & HT_MASK;
    while (g_ht_used[h]) {
        if (g_ht[h].r == r && g_ht[h].m == m &&
            g_ht[h].unused == unused && g_ht[h].remaining == remaining)
            return 1;
        h = (h + 1) & HT_MASK;
    }
    return 0;
}

static void ht_insert(i64 r, i64 m, i64 unused, i64 remaining) {
    u64 h = state_hash(r, m, unused, remaining) & HT_MASK;
    while (g_ht_used[h]) {
        if (g_ht[h].r == r && g_ht[h].m == m &&
            g_ht[h].unused == unused && g_ht[h].remaining == remaining)
            return;
        h = (h + 1) & HT_MASK;
    }
    g_ht[h].r = r; g_ht[h].m = m;
    g_ht[h].unused = unused; g_ht[h].remaining = remaining;
    g_ht_used[h] = 1;
}

/* ---------- DFS ---------- */

static i64 g_M;
static i64 g_residues[60000];
static int g_num_residues;

static int candidates_offsets(int Lval, i64 unused_mask, int target, int step, int *out) {
    int n = 0;
    int j = target;
    while (j < Lval) {
        if ((unused_mask >> j) & 1)
            out[n++] = j;
        j += step;
    }
    return n;
}

static int pick_index_mrv(int Lval, i64 r, i64 m, i64 unused_mask, i64 remaining_mask,
                          int *best_i_out, int *best_cands_out) {
    int best_i = 0;
    int best_count = 1000000000;
    int best_cands[Lval];
    int best_n = 0;

    for (int i = Lval; i >= 1; i--) {
        if (!((remaining_mask >> (i - 1)) & 1))
            continue;
        i64 g = gcd64(m, i);
        i64 target = (((-r) % g) + g) % g;
        int cands[Lval];
        int n = candidates_offsets(Lval, unused_mask, (int)target, (int)g, cands);
        if (n == 0)
            return 0;
        if (n < best_count || (n == best_count && i > best_i)) {
            best_count = n;
            best_i = i;
            memcpy(best_cands, cands, n * sizeof(int));
            best_n = n;
        }
    }

    *best_i_out = best_i;
    memcpy(best_cands_out, best_cands, best_n * sizeof(int));
    return best_n;
}

static void dfs(i64 r, i64 m, i64 unused_mask, i64 remaining_mask) {
    r %= m;
    if (ht_contains(r, m, unused_mask, remaining_mask))
        return;
    ht_insert(r, m, unused_mask, remaining_mask);

    if (remaining_mask == 0) {
        g_residues[g_num_residues++] = r;
        return;
    }

    int best_i, best_cands[L];
    int best_n = pick_index_mrv(L, r, m, unused_mask, remaining_mask, &best_i, best_cands);
    if (best_n == 0)
        return;

    i64 remaining2 = remaining_mask & ~((i64)1 << (best_i - 1));

    for (int idx = 0; idx < best_n; idx++) {
        int j = best_cands[idx];
        i64 r2_val = (((-j) % best_i) + best_i) % best_i;
        i64 new_r, new_m;
        if (crt_merge(r, m, r2_val, best_i, &new_r, &new_m)) {
            dfs(new_r, new_m, unused_mask & ~((i64)1 << j), remaining2);
        }
    }
}

/* ---------- Sorting ---------- */

static int cmp_i64(const void *a, const void *b) {
    i64 va = *(const i64 *)a;
    i64 vb = *(const i64 *)b;
    if (va < vb) return -1;
    if (va > vb) return 1;
    return 0;
}

/* ---------- Main ---------- */

long long p896_native(void) {
    g_M = lcm_upto(L);

    i64 all_offsets = ((i64)1 << L) - 1;
    i64 all_indices = ((i64)1 << L) - 1;

    g_num_residues = 0;
    dfs(0, 1, all_offsets, all_indices);

    /* Sort and deduplicate */
    qsort(g_residues, g_num_residues, sizeof(i64), cmp_i64);
    int j = 0;
    for (int i = 0; i < g_num_residues; i++) {
        if (i == 0 || g_residues[i] != g_residues[i - 1])
            g_residues[j++] = g_residues[i];
    }
    g_num_residues = j;

    /* Convert to starts: 0 -> M */
    for (int i = 0; i < g_num_residues; i++) {
        if (g_residues[i] == 0)
            g_residues[i] = g_M;
    }

    /* Sort again */
    qsort(g_residues, g_num_residues, sizeof(i64), cmp_i64);

    return g_residues[36 - 1];
}
