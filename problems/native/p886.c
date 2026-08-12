// Project Euler 886: Coprime Permutations.
// For n=2m, counts Hamiltonian cycles in a balanced bipartite coprimality graph
// using determinant/permanent inclusion-exclusion with grouped classes.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long i64;
typedef unsigned long long u64;
typedef uint8_t u8;
typedef __int128 i128;

#define MOD 83456729LL
#define MAXN 34

static i64 mulmod(i64 a, i64 b, i64 m) {
    return (i64)((i128)a * b % m);
}

static i64 powmod(i64 a, i64 e, i64 m) {
    i64 r = 1 % m; a %= m; if (a < 0) a += m;
    while (e > 0) {
        if (e & 1) r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}

static i64 gcd_ll(i64 a, i64 b) {
    while (b) { i64 t = a % b; a = b; b = t; }
    return a;
}

// Binomial table
static i64 comb[MAXN + 1][MAXN + 1];

static void init_comb(int n) {
    for (int i = 0; i <= n; i++) {
        comb[i][0] = comb[i][i] = 1;
        for (int j = 1; j < i; j++)
            comb[i][j] = comb[i-1][j-1] + comb[i-1][j];
    }
}

// Determinant mod using Gaussian elimination
static i64 determinant_mod(i64 *mat, int size, i64 mod) {
    if (size == 0) return 1;
    i64 a[MAXN][MAXN];
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++) {
            a[i][j] = mat[i * size + j] % mod;
            if (a[i][j] < 0) a[i][j] += mod;
        }
    i64 det = 1;
    for (int col = 0; col < size; col++) {
        int pivot = -1;
        for (int row = col; row < size; row++) {
            if (a[row][col] % mod != 0) { pivot = row; break; }
        }
        if (pivot < 0) return 0;
        if (pivot != col) {
            for (int j = 0; j < size; j++) {
                i64 t = a[col][j]; a[col][j] = a[pivot][j]; a[pivot][j] = t;
            }
            det = -det;
        }
        i64 pv = a[col][col] % mod;
        if (pv < 0) pv += mod;
        det = mulmod(det, pv, mod);
        i64 inv = powmod(pv, mod - 2, mod);
        for (int row = col + 1; row < size; row++) {
            if (a[row][col] % mod != 0) {
                i64 factor = mulmod(a[row][col], inv, mod);
                for (int j = col; j < size; j++) {
                    a[row][j] = (a[row][j] - mulmod(factor, a[col][j], mod)) % mod;
                    if (a[row][j] < 0) a[row][j] += mod;
                }
            }
        }
    }
    det %= mod;
    if (det < 0) det += mod;
    return det;
}

// ---- Global context for grouped_permanent ----
static int g_row_class_count, g_col_class_count;
static int g_row_class_rep[MAXN], g_row_class_count_arr[MAXN];
static int g_col_class_rep[MAXN], g_col_class_count_arr[MAXN];
static int g_class_entry[MAXN][MAXN];  // [row_class][col_class]
static int g_row_supports[MAXN];       // support bitmask over col classes

// permanent cache: key = (row_selected << 20) | col_selected
#define PERM_CACHE_CAP (1 << 20)
#define PERM_CACHE_MASK (PERM_CACHE_CAP - 1)
typedef struct { u64 key; i64 val; u8 used; } PermEntry;
static PermEntry g_perm_cache[PERM_CACHE_CAP];

static u64 perm_hash(u64 key) {
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    return key & PERM_CACHE_MASK;
}

static int perm_cache_lookup(u64 key, i64 *out) {
    u64 h = perm_hash(key);
    while (g_perm_cache[h].used) {
        if (g_perm_cache[h].key == key) { *out = g_perm_cache[h].val; return 1; }
        h = (h + 1) & PERM_CACHE_MASK;
    }
    return 0;
}

static void perm_cache_insert(u64 key, i64 val) {
    u64 h = perm_hash(key);
    while (g_perm_cache[h].used) {
        if (g_perm_cache[h].key == key) { g_perm_cache[h].val = val; return; }
        h = (h + 1) & PERM_CACHE_MASK;
    }
    g_perm_cache[h].used = 1;
    g_perm_cache[h].key = key;
    g_perm_cache[h].val = val;
}

// Recursive context for grouped_permanent
static int g_active_col_classes[MAXN];
static int g_active_col_counts[MAXN];
static int g_active_count;
static int g_supports[256];   // support bitmasks
static int g_row_counts_arr[256];
static int g_row_sums[256];
static int g_num_supports;
static int g_row_total;
static i64 g_total;

static void rec_perm(int col_pos, int chosen, i64 coeff) {
    if (col_pos == g_active_count) {
        i64 product = coeff;
        for (int i = 0; i < g_num_supports; i++) {
            if (g_row_sums[i] == 0) return;
            product = mulmod(product, powmod(g_row_sums[i], g_row_counts_arr[i], MOD), MOD);
        }
        if (chosen & 1) g_total -= product;
        else g_total += product;
        return;
    }
    int class_size = g_active_col_counts[col_pos];
    // Find affected indices
    int affected[256];
    int n_affected = 0;
    for (int idx = 0; idx < g_num_supports; idx++) {
        if ((g_supports[idx] >> col_pos) & 1)
            affected[n_affected++] = idx;
    }
    for (int take = 0; take <= class_size; take++) {
        if (take) {
            for (int i = 0; i < n_affected; i++)
                g_row_sums[affected[i]] += take;
        }
        rec_perm(col_pos + 1, chosen + take,
                 mulmod(coeff, comb[class_size][take] % MOD, MOD));
        if (take) {
            for (int i = 0; i < n_affected; i++)
                g_row_sums[affected[i]] -= take;
        }
    }
}

static i64 grouped_permanent(int row_selected, int col_selected) {
    u64 key = ((u64)row_selected << 20) | (u64)col_selected;
    i64 cached;
    if (perm_cache_lookup(key, &cached)) return cached;

    // Build active col classes
    g_active_count = 0;
    for (int cc = 0; cc < g_col_class_count; cc++) {
        int remaining = g_col_class_count_arr[cc] - ((col_selected >> cc) & 1);
        if (remaining > 0) {
            g_active_col_classes[g_active_count] = cc;
            g_active_col_counts[g_active_count] = remaining;
            g_active_count++;
        }
    }

    // Build row groups using a simple array indexed by support bitmask
    // The distinguished row (row 0) always has full support and count 1.
    int row_group_counts[256];
    memset(row_group_counts, 0, sizeof(row_group_counts));
    row_group_counts[(1 << g_active_count) - 1] = 1;
    for (int rc = 0; rc < g_row_class_count; rc++) {
        int remaining = g_row_class_count_arr[rc] - ((row_selected >> rc) & 1);
        if (remaining <= 0) continue;
        int support = 0;
        int orig = g_row_supports[rc];
        for (int nc = 0; nc < g_active_count; nc++) {
            int old_col = g_active_col_classes[nc];
            if ((orig >> old_col) & 1)
                support |= 1 << nc;
        }
        row_group_counts[support] += remaining;
    }

    // Collect supports and counts
    g_num_supports = 0;
    g_row_total = 0;
    for (int s = 0; s < 256; s++) {
        if (row_group_counts[s] > 0) {
            g_supports[g_num_supports] = s;
            g_row_counts_arr[g_num_supports] = row_group_counts[s];
            g_row_sums[g_num_supports] = 0;
            g_row_total += row_group_counts[s];
            g_num_supports++;
        }
    }

    g_total = 0;
    rec_perm(0, 0, 1);
    if (g_row_total & 1) g_total = -g_total;
    g_total %= MOD;
    if (g_total < 0) g_total += MOD;

    perm_cache_insert(key, g_total);
    return g_total;
}

static i64 P_func(int n, i64 mod) {
    // Clear permanent cache from any previous call
    memset(g_perm_cache, 0, sizeof(g_perm_cache));
    if (n < 2) return 0;
    int m = n / 2;

    // Build coprimality matrix
    int odds[MAXN], evens[MAXN];
    for (int i = 0; i < m; i++) odds[i] = 2 * i + 1;
    for (int i = 0; i < m; i++) evens[i] = 2 * (i + 1);

    int matrix[MAXN][MAXN];
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            matrix[i][j] = (gcd_ll(odds[i], evens[j]) == 1) ? 1 : 0;

    // Row masks for rows 1..m-1
    int row_masks[MAXN];
    for (int row = 1; row < m; row++) {
        int mask = 0;
        for (int col = 0; col < m; col++)
            if (matrix[row][col]) mask |= 1 << col;
        row_masks[row - 1] = mask;
    }

    // Col masks for all columns
    int col_masks[MAXN];
    for (int col = 0; col < m; col++) {
        int mask = 0;
        for (int row = 0; row < m; row++)
            if (matrix[row][col]) mask |= 1 << row;
        col_masks[col] = mask;
    }

    // Group rows by mask (preserve first-seen order)
    g_row_class_count = 0;
    int *row_seen = (int *)malloc(131072 * sizeof(int));
    memset(row_seen, -1, 131072 * sizeof(int));
    for (int idx = 0; idx < m - 1; idx++) {
        int mask = row_masks[idx];
        if (row_seen[mask] < 0) {
            row_seen[mask] = g_row_class_count;
            g_row_class_rep[g_row_class_count] = idx;
            g_row_class_count_arr[g_row_class_count] = 1;
            g_row_class_count++;
        } else {
            g_row_class_count_arr[row_seen[mask]]++;
        }
    }

    // Group cols by mask
    g_col_class_count = 0;
    int *col_seen = (int *)malloc(131072 * sizeof(int));
    memset(col_seen, -1, 131072 * sizeof(int));
    for (int idx = 0; idx < m; idx++) {
        int mask = col_masks[idx];
        if (col_seen[mask] < 0) {
            col_seen[mask] = g_col_class_count;
            g_col_class_rep[g_col_class_count] = idx;
            g_col_class_count_arr[g_col_class_count] = 1;
            g_col_class_count++;
        } else {
            g_col_class_count_arr[col_seen[mask]]++;
        }
    }

    // class_entry: [row_class][col_class]
    for (int rc = 0; rc < g_row_class_count; rc++)
        for (int cc = 0; cc < g_col_class_count; cc++)
            g_class_entry[rc][cc] = matrix[g_row_class_rep[rc] + 1][g_col_class_rep[cc]];

    // row_supports: support over col classes
    for (int rc = 0; rc < g_row_class_count; rc++) {
        int support = 0;
        for (int cc = 0; cc < g_col_class_count; cc++) {
            if (matrix[g_row_class_rep[rc] + 1][g_col_class_rep[cc]])
                support |= 1 << cc;
        }
        g_row_supports[rc] = support;
    }

    init_comb(m);

    // Enumerate combinations
    // For each size, enumerate all subsets of row classes and col classes
    i64 answer = 0;
    int rcc = g_row_class_count, ccc = g_col_class_count;

    // Generate all subsets of row classes and col classes
    // For each subset size, pair them up
    for (int size = 0; size <= ccc && size <= rcc; size++) {
        i64 sign = (size & 1) ? -1 : 1;

        // Enumerate row subsets of given size
        int row_subset[20];
        // Use combination enumeration
        // Initialize
        for (int i = 0; i < size; i++) row_subset[i] = i;
        while (1) {
            // Process current row subset
            int row_mask = 0;
            i64 row_mult = 1;
            for (int i = 0; i < size; i++) {
                row_mask |= 1 << row_subset[i];
                row_mult = mulmod(row_mult, g_row_class_count_arr[row_subset[i]], mod);
            }

            // Enumerate col subsets of same size
            int col_subset[20];
            for (int i = 0; i < size; i++) col_subset[i] = i;
            while (1) {
                int col_mask = 0;
                i64 col_mult = 1;
                for (int i = 0; i < size; i++) {
                    col_mask |= 1 << col_subset[i];
                    col_mult = mulmod(col_mult, g_col_class_count_arr[col_subset[i]], mod);
                }

                // Build det matrix
                i64 det_mat[400];  // 20x20 max
                for (int i = 0; i < size; i++)
                    for (int j = 0; j < size; j++)
                        det_mat[i * size + j] = g_class_entry[row_subset[i]][col_subset[j]];

                i64 det = determinant_mod(det_mat, size, mod);
                if (det != 0) {
                    i64 perm = grouped_permanent(row_mask, col_mask);
                    i64 term = mulmod(row_mult, col_mult, mod);
                    term = mulmod(term, mulmod(det, det, mod), mod);
                    term = mulmod(term, mulmod(perm, perm, mod), mod);
                    answer += sign * term;
                }

                // Next col subset
                if (size == 0) break;
                int pos = size - 1;
                while (pos >= 0 && col_subset[pos] == ccc - size + pos) pos--;
                if (pos < 0) break;
                col_subset[pos]++;
                for (int i = pos + 1; i < size; i++)
                    col_subset[i] = col_subset[i-1] + 1;
            }

            // Next row subset
            if (size == 0) break;
            int pos = size - 1;
            while (pos >= 0 && row_subset[pos] == rcc - size + pos) pos--;
            if (pos < 0) break;
            row_subset[pos]++;
            for (int i = pos + 1; i < size; i++)
                row_subset[i] = row_subset[i-1] + 1;
        }
    }

    answer %= mod;
    if (answer < 0) answer += mod;
    free(row_seen);
    free(col_seen);
    return answer;
}

long long p886_native(void) {
    return (long long)P_func(34, MOD);
}
