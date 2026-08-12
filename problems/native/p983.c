// Project Euler 983: Parity family of circles on lattice points.
// Find minimum radius^2 m for k=10 vectors giving 512 circles, 512 harmony points.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Encode (x,y) as int64 for hash set/map
static long long encode(int x, int y) {
    return ((long long)(x + 100000) << 32) | (unsigned int)(y + 100000);
}

// Simple hash set for long long keys
typedef struct {
    long long *keys;
    int *used;  // 0=empty, 1=occupied
    int cap;
    int count;
} HashSet;

static void hs_init(HashSet *hs, int cap) {
    hs->cap = cap;
    hs->count = 0;
    hs->keys = calloc(cap, sizeof(long long));
    hs->used = calloc(cap, sizeof(int));
}

static void hs_free(HashSet *hs) {
    free(hs->keys);
    free(hs->used);
}

static unsigned int hash_ll(long long key, int cap) {
    return (unsigned int)(key % cap + cap) % cap;
}

static void hs_insert(HashSet *hs, long long key) {
    if (hs->count * 2 > hs->cap) {
        // resize
        int old_cap = hs->cap;
        long long *old_keys = hs->keys;
        int *old_used = hs->used;
        hs->cap = old_cap * 2;
        hs->keys = calloc(hs->cap, sizeof(long long));
        hs->used = calloc(hs->cap, sizeof(int));
        hs->count = 0;
        for (int i = 0; i < old_cap; i++) {
            if (old_used[i]) {
                unsigned int h = hash_ll(old_keys[i], hs->cap);
                while (hs->used[h]) h = (h + 1) % hs->cap;
                hs->keys[h] = old_keys[i];
                hs->used[h] = 1;
                hs->count++;
            }
        }
        free(old_keys);
        free(old_used);
    }
    unsigned int h = hash_ll(key, hs->cap);
    while (hs->used[h]) {
        if (hs->keys[h] == key) return;
        h = (h + 1) % hs->cap;
    }
    hs->keys[h] = key;
    hs->used[h] = 1;
    hs->count++;
}

static int hs_contains(HashSet *hs, long long key) {
    unsigned int h = hash_ll(key, hs->cap);
    while (hs->used[h]) {
        if (hs->keys[h] == key) return 1;
        h = (h + 1) % hs->cap;
    }
    return 0;
}

// Simple hash map for long long -> int
typedef struct {
    long long *keys;
    int *vals;
    int *used;
    int cap;
} HashMap;

static void hm_init(HashMap *hm, int cap) {
    hm->cap = cap;
    hm->keys = calloc(cap, sizeof(long long));
    hm->vals = calloc(cap, sizeof(int));
    hm->used = calloc(cap, sizeof(int));
}

static void hm_free(HashMap *hm) {
    free(hm->keys);
    free(hm->vals);
    free(hm->used);
}

static int hm_get(HashMap *hm, long long key, int *found) {
    unsigned int h = hash_ll(key, hm->cap);
    while (hm->used[h]) {
        if (hm->keys[h] == key) { *found = 1; return hm->vals[h]; }
        h = (h + 1) % hm->cap;
    }
    *found = 0;
    return 0;
}

static void hm_set(HashMap *hm, long long key, int val) {
    unsigned int h = hash_ll(key, hm->cap);
    while (hm->used[h]) {
        if (hm->keys[h] == key) { hm->vals[h] = val; return; }
        h = (h + 1) % hm->cap;
    }
    hm->keys[h] = key;
    hm->vals[h] = val;
    hm->used[h] = 1;
}

static int isqrt_int(int n) {
    if (n < 0) return 0;
    int r = (int)sqrt((double)n);
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r;
}

static int antipodal_pair_count(int m) {
    int x = m;
    while (x % 2 == 0) x /= 2;
    int product = 1;
    for (int p = 3; p * p <= x; p += 2) {
        if (x % p == 0) {
            int exp = 0;
            while (x % p == 0) { x /= p; exp++; }
            if (p % 4 == 1) product *= exp + 1;
            else if (exp & 1) return 0;
        }
    }
    if (x > 1) {
        if (x % 4 == 1) product *= 2;
        else if (x % 4 == 3) return 0;
    }
    return 2 * product;
}

#define MAX_POINTS 100

static int lattice_points_on_circle(int m, int px[], int py[]) {
    int lim = isqrt_int(m);
    int cnt = 0;
    for (int x = -lim; x <= lim; x++) {
        int y2 = m - x * x;
        if (y2 < 0) continue;
        int y = isqrt_int(y2);
        if (y * y == y2) {
            px[cnt] = x; py[cnt] = y; cnt++;
            if (y) { px[cnt] = x; py[cnt] = -y; cnt++; }
        }
    }
    return cnt;
}

static int cmp_points(const void *a, const void *b) {
    const int *pa = (const int *)a;
    const int *pb = (const int *)b;
    if (pa[0] != pb[0]) return pa[0] - pb[0];
    return pa[1] - pb[1];
}

// Global state for the search
static int circle_px[MAX_POINTS], circle_py[MAX_POINTS];
static int num_circle_points;
static HashSet deltas;
static int num_pairs;
static int pair_px[MAX_POINTS], pair_py[MAX_POINTS];  // first element of each pair
static int pair_nx[MAX_POINTS], pair_ny[MAX_POINTS];  // negation

static int has_unit_coordinate(int npts, int px[], int py[]) {
    for (int i = 0; i < npts; i++) {
        if (abs(px[i]) == 1 || abs(py[i]) == 1) return 1;
    }
    return 0;
}

static void build_displacement_set() {
    hs_init(&deltas, 1024);
    for (int i = 0; i < num_circle_points; i++) {
        for (int j = 0; j < num_circle_points; j++) {
            if (i != j) {
                hs_insert(&deltas, encode(circle_px[i] - circle_px[j],
                                          circle_py[i] - circle_py[j]));
            }
        }
    }
}

static int passes_four_vector_prune(int sel_x[], int sel_y[], int sel_len,
                                     int vx, int vy) {
    if (sel_len < 3) return 1;
    for (int i = 0; i < sel_len - 2; i++) {
        int ax = sel_x[i], ay = sel_y[i];
        for (int j = i + 1; j < sel_len - 1; j++) {
            int bx = sel_x[j], by = sel_y[j];
            for (int k = j + 1; k < sel_len; k++) {
                int cx = sel_x[k], cy = sel_y[k];
                for (int sa = 1; sa >= -1; sa -= 2) {
                    int x1 = vx + sa * ax;
                    int y1 = vy + sa * ay;
                    for (int sb = 1; sb >= -1; sb -= 2) {
                        int x2 = x1 + sb * bx;
                        int y2 = y1 + sb * by;
                        // +cx
                        int x = x2 + cx, y = y2 + cy;
                        if ((x || y) && hs_contains(&deltas, encode(x, y))) return 0;
                        // -cx
                        x = x2 - cx; y = y2 - cy;
                        if ((x || y) && hs_contains(&deltas, encode(x, y))) return 0;
                    }
                }
            }
        }
    }
    return 1;
}

// Even masks for k bits
static int even_masks[1024];
static int num_even_masks;

static void compute_even_masks(int k) {
    num_even_masks = 0;
    for (int mask = 0; mask < (1 << k); mask++) {
        if (__builtin_popcount(mask) % 2 == 0) {
            even_masks[num_even_masks++] = mask;
        }
    }
}

static void centers_from_vectors(int vx[], int vy[], int k, int masks[], int nmasks,
                                  int cx[], int cy[]) {
    for (int mi = 0; mi < nmasks; mi++) {
        int mask = masks[mi];
        int x = 0, y = 0;
        for (int i = 0; i < k; i++) {
            if (mask & (1 << i)) { x += vx[i]; y += vy[i]; }
        }
        cx[mi] = x; cy[mi] = y;
    }
}

static int quick_harmony_count_equals_n(int ncen, int cx[], int cy[], int n) {
    HashMap counts;
    hm_init(&counts, 16384);
    int harmony = 0;
    for (int ci = 0; ci < ncen; ci++) {
        for (int vi = 0; vi < num_circle_points; vi++) {
            long long key = encode(cx[ci] + circle_px[vi], cy[ci] + circle_py[vi]);
            int found;
            int cur = hm_get(&counts, key, &found);
            if (!found) {
                hm_set(&counts, key, 1);
            } else if (cur == 1) {
                hm_set(&counts, key, 2);
                harmony++;
                if (harmony > n) { hm_free(&counts); return 0; }
            } else {
                hm_set(&counts, key, cur + 1);
            }
        }
    }
    hm_free(&counts);
    return harmony == n;
}

// Union-Find
static int uf_parent[1024], uf_size[1024];

static int uf_find(int x) {
    while (uf_parent[x] != x) {
        uf_parent[x] = uf_parent[uf_parent[x]];
        x = uf_parent[x];
    }
    return x;
}

static void uf_union(int a, int b) {
    int ra = uf_find(a), rb = uf_find(b);
    if (ra == rb) return;
    if (uf_size[ra] < uf_size[rb]) { int t = ra; ra = rb; rb = t; }
    uf_parent[rb] = ra;
    uf_size[ra] += uf_size[rb];
}

static int strict_perfect_check(int ncen, int cx[], int cy[], int n) {
    // Check no tangent circle pairs
    HashSet center_set;
    hs_init(&center_set, 1024);
    for (int i = 0; i < ncen; i++)
        hs_insert(&center_set, encode(cx[i], cy[i]));

    for (int i = 0; i < ncen; i++) {
        for (int vi = 0; vi < num_circle_points; vi++) {
            long long other = encode(cx[i] + 2 * circle_px[vi],
                                      cy[i] + 2 * circle_py[vi]);
            if (hs_contains(&center_set, other)) {
                long long ci_enc = encode(cx[i], cy[i]);
                if (ci_enc < other) {
                    hs_free(&center_set);
                    return 0;
                }
            }
        }
    }

    // Build harmony-point incidences
    HashMap point_to_centers;
    hm_init(&point_to_centers, 16384);
    // We need a list of centers per point. Use a different approach:
    // Store comma-separated center indices as val (packed).
    // Actually, let's use a simpler approach: for each point, store first center index,
    // and a linked list.

    // Let's use a hash map from point key to a dynamic list of center indices.
    // Since we can't easily do that, let's use a different approach:
    // For each pair of circles, check if they share a point.

    // Actually, let me just count harmony points and build components.
    // Reuse the hash map to store count, then iterate again for components.

    // Count points
    HashMap pcounts;
    hm_init(&pcounts, 16384);

    for (int ci = 0; ci < ncen; ci++) {
        for (int vi = 0; vi < num_circle_points; vi++) {
            long long key = encode(cx[ci] + circle_px[vi], cy[ci] + circle_py[vi]);
            int found;
            int cur = hm_get(&pcounts, key, &found);
            if (!found) hm_set(&pcounts, key, 1);
            else hm_set(&pcounts, key, cur + 1);
        }
    }

    // Count harmony points (count >= 2)
    int harmony_count = 0;
    // We need to identify which points are harmony points and which circles they connect.
    // Let's do it differently: for each pair of circles, find shared points.

    // Actually, let's just check all pairs of circles for shared lattice points.
    // With 512 circles and 24 points each, that's 512*512/2 = 131072 pairs,
    // each checking 24*24 = 576 point pairs. That's 75M operations. Might be slow.

    // Better: build a hash map from point -> list of center indices.
    // Let me use a hash map with chaining.

    // Let me use a simpler approach: create a list of (point_key, center_idx) pairs,
    // sort by point_key, then group.

    int total_entries = ncen * num_circle_points;
    long long *entry_keys = malloc(total_entries * sizeof(long long));
    int *entry_centers = malloc(total_entries * sizeof(int));

    int idx = 0;
    for (int ci = 0; ci < ncen; ci++) {
        for (int vi = 0; vi < num_circle_points; vi++) {
            entry_keys[idx] = encode(cx[ci] + circle_px[vi], cy[ci] + circle_py[vi]);
            entry_centers[idx] = ci;
            idx++;
        }
    }

    // Sort by key (simple insertion sort or qsort with paired data)
    // Use qsort on indices
    // Actually, let's sort using a simple approach: create index array, sort that
    int *sort_idx = malloc(total_entries * sizeof(int));
    for (int i = 0; i < total_entries; i++) sort_idx[i] = i;

    // Simple sort: bubble sort is too slow. Use qsort with global arrays.
    // Let's use a different approach: since we have the hash map pcounts,
    // let's just iterate through all entries and for each harmony point,
    // find all centers that have it.

    // Actually, let me just do the O(n^2 * p^2) approach but optimized.
    // For each pair of circles, check if they share any lattice point.
    // A shared point means: cx[i]+vx == cx[j]+vx' and cy[i]+vy == cy[j]+vy'
    // => cx[i]-cx[j] == vx'-vx and cy[i]-cy[j] == vy'-vy
    // So the difference of centers must be a displacement between circle points.

    // We already have the displacement set! So two circles share a point iff
    // their center difference is in the displacement set.

    int harmony = 0;
    for (int i = 0; i < ncen; i++) {
        for (int j = i + 1; j < ncen; j++) {
            int dx = cx[i] - cx[j];
            int dy = cy[i] - cy[j];
            if (hs_contains(&deltas, encode(dx, dy))) {
                // They share at least one point
                // But we need to count the number of shared points, not just whether they share
                // Actually, the harmony count is the number of points touched by >= 2 circles,
                // not the number of pairs. Let me reconsider.

                // Actually, let me go back to the sorted approach.
            }
        }
    }

    // Let me use the sorted entries approach properly.
    // Sort entry_keys and entry_centers together using a simple merge sort or
    // just use the hash map approach differently.

    // Let me use a hash map that stores a linked list of center indices.
    // I'll use the hash map to store the head of a linked list, and
    // arrays for the linked list.

    int *link_head = malloc(total_entries * sizeof(int));  // index in link_next
    int *link_next = malloc(total_entries * sizeof(int));
    int *link_center = malloc(total_entries * sizeof(int));
    long long *link_key = malloc(total_entries * sizeof(long long));

    HashMap link_map;  // key -> index in link arrays
    hm_init(&link_map, 16384);

    int link_count = 0;
    for (int ci = 0; ci < ncen; ci++) {
        for (int vi = 0; vi < num_circle_points; vi++) {
            long long key = encode(cx[ci] + circle_px[vi], cy[ci] + circle_py[vi]);
            int found;
            int head = hm_get(&link_map, key, &found);
            if (!found) {
                link_key[link_count] = key;
                link_center[link_count] = ci;
                link_next[link_count] = -1;
                hm_set(&link_map, key, link_count);
                link_count++;
            } else {
                link_key[link_count] = key;
                link_center[link_count] = ci;
                link_next[link_count] = head;
                hm_set(&link_map, key, link_count);
                link_count++;
            }
        }
    }

    // Now count harmony points and build union-find
    for (int i = 0; i < n; i++) { uf_parent[i] = i; uf_size[i] = 1; }

    harmony = 0;
    // Iterate through all unique keys in the hash map
    // Since we can't easily iterate a hash map, let's iterate through link entries
    // and group by key.

    // Actually, let's iterate through the link_map's occupied entries.
    for (int h = 0; h < link_map.cap; h++) {
        if (!link_map.used[h]) continue;
        long long key = link_map.keys[h];
        int head = link_map.vals[h];

        // Follow the linked list to count centers
        int count = 0;
        int first_center = -1;
        int idx2 = head;
        while (idx2 != -1) {
            count++;
            if (first_center == -1) first_center = link_center[idx2];
            idx2 = link_next[idx2];
        }

        if (count >= 2) {
            harmony++;
            // Union all centers sharing this point
            idx2 = head;
            int base = link_center[idx2];
            idx2 = link_next[idx2];
            while (idx2 != -1) {
                uf_union(base, link_center[idx2]);
                idx2 = link_next[idx2];
            }
        }
    }

    if (harmony != n) {
        hm_free(&link_map);
        free(link_head); free(link_next); free(link_center); free(link_key);
        free(entry_keys); free(entry_centers); free(sort_idx);
        hm_free(&pcounts);
        hs_free(&center_set);
        return 0;
    }

    // Check all circles are connected
    int root = uf_find(0);
    int connected = 1;
    for (int i = 1; i < n; i++) {
        if (uf_find(i) != root) { connected = 0; break; }
    }

    hm_free(&link_map);
    free(link_head); free(link_next); free(link_center); free(link_key);
    free(entry_keys); free(entry_centers); free(sort_idx);
    hm_free(&pcounts);
    hs_free(&center_set);

    return connected;
}

// DFS state
static int sel_x[16], sel_y[16];
static int dfs_k, dfs_n;
static int dfs_masks[1024], dfs_nmasks;

static int dfs(int start, int sel_len) {
    if (sel_len == dfs_k) {
        int cx[1024], cy[1024];
        centers_from_vectors(sel_x, sel_y, dfs_k, dfs_masks, dfs_nmasks, cx, cy);
        if (!quick_harmony_count_equals_n(dfs_nmasks, cx, cy, dfs_n)) return 0;
        return strict_perfect_check(dfs_nmasks, cx, cy, dfs_n);
    }

    int needed = dfs_k - sel_len;
    for (int pi = start; pi < num_pairs - needed + 1; pi++) {
        // First pair: only one orientation; others: both
        int n_choices = (sel_len == 0) ? 1 : 2;
        for (int ci = 0; ci < n_choices; ci++) {
            int vx, vy;
            if (ci == 0) { vx = pair_px[pi]; vy = pair_py[pi]; }
            else { vx = pair_nx[pi]; vy = pair_ny[pi]; }

            if (!passes_four_vector_prune(sel_x, sel_y, sel_len, vx, vy))
                continue;

            sel_x[sel_len] = vx;
            sel_y[sel_len] = vy;
            if (dfs(pi + 1, sel_len + 1)) return 1;
        }
    }
    return 0;
}

static int has_valid_oriented_vectors(int k, int n) {
    dfs_k = k;
    dfs_n = n;
    compute_even_masks(k);
    dfs_nmasks = num_even_masks;
    for (int i = 0; i < num_even_masks; i++) dfs_masks[i] = even_masks[i];
    return dfs(0, 0);
}

static int find_min_radius_sq(int k, int m_limit, int filtered) {
    int n = 1 << (k - 1);

    for (int m = 1; m <= m_limit; m++) {
        int p = antipodal_pair_count(m);
        if (p < k) continue;
        if (filtered && p != k && p != k + 2) continue;

        num_circle_points = lattice_points_on_circle(m, circle_px, circle_py);

        if (filtered) {
            if (!has_unit_coordinate(num_circle_points, circle_px, circle_py))
                continue;
        }

        // Build opposite pairs
        // Sort points
        int sorted_pts[MAX_POINTS][2];
        for (int i = 0; i < num_circle_points; i++) {
            sorted_pts[i][0] = circle_px[i];
            sorted_pts[i][1] = circle_py[i];
        }
        qsort(sorted_pts, num_circle_points, sizeof(int) * 2, cmp_points);

        // Group into antipodal pairs
        char used[MAX_POINTS] = {0};
        num_pairs = 0;
        for (int i = 0; i < num_circle_points; i++) {
            if (used[i]) continue;
            int nx = -sorted_pts[i][0], ny = -sorted_pts[i][1];
            // Find negation in sorted array
            int lo = 0, hi = num_circle_points - 1, found = -1;
            while (lo <= hi) {
                int mid = (lo + hi) / 2;
                if (sorted_pts[mid][0] < nx || (sorted_pts[mid][0] == nx && sorted_pts[mid][1] < ny))
                    lo = mid + 1;
                else if (sorted_pts[mid][0] > nx || (sorted_pts[mid][0] == nx && sorted_pts[mid][1] > ny))
                    hi = mid - 1;
                else { found = mid; break; }
            }
            if (found == -1) continue;  // shouldn't happen
            used[i] = 1;
            used[found] = 1;
            pair_px[num_pairs] = sorted_pts[i][0];
            pair_py[num_pairs] = sorted_pts[i][1];
            pair_nx[num_pairs] = nx;
            pair_ny[num_pairs] = ny;
            num_pairs++;
        }

        if (num_pairs != p) continue;

        build_displacement_set();

        if (has_valid_oriented_vectors(k, n)) {
            hs_free(&deltas);
            return m;
        }
        hs_free(&deltas);
    }
    return -1;  // not found
}

long long p983_native(void) {
    return find_min_radius_sq(10, 20000, 1);
}
