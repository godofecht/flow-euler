/* Project Euler 933: Paper Cutting.
 *
 * Two-player impartial game on integer rectangles. A move cuts one
 * rectangle once vertically and once horizontally into four smaller
 * rectangles (all sides positive integers). No move means lose.
 *
 * Let C(w,h) be the number of winning first moves from a w x h rectangle.
 * Let D(W,H) = sum_{w=2..W} sum_{h=2..H} C(w,h).
 *
 * We compute D(123, 1234567).
 *
 * Approach:
 * 1) Compute Sprague-Grundy numbers G(w,h) (nimbers).
 * 2) For each fixed width w, prove G(w,h) becomes constant beyond some
 *    index s by detecting it is constant on [s, 2s].
 * 3) Count D(W,H) via quadruple counting over ordered splits.
 *
 * Ported from the Python reference solver. The answer fits in int64.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ----------------- GrundyComputer state ----------------- */

static int *B_arr;        /* B[w] = computed height limit for width w */
static int *g_const_arr;  /* g_const[w] = stabilization start index */
static int *const_arr;    /* const[w] = constant nimber after stabilization */
static int **G_arr;       /* G_arr[w][h] = nimber for w x h, h=0..B[w] */
static int *visited;      /* stamp-based visited array for mex */
static int vis_size;
static int stamp;
static int prev_maxg;
static int current_w;
static int W_alloc;

static void grundy_init(void) {
    W_alloc = 4;
    B_arr = calloc(W_alloc, sizeof(int));
    g_const_arr = calloc(W_alloc, sizeof(int));
    const_arr = calloc(W_alloc, sizeof(int));
    G_arr = calloc(W_alloc, sizeof(int *));
    B_arr[1] = 1;
    g_const_arr[1] = 1;
    G_arr[1] = calloc(2, sizeof(int));
    G_arr[1][0] = 0;
    G_arr[1][1] = 0;
    vis_size = 32768;
    visited = calloc(vis_size, sizeof(int));
    stamp = 1;
    prev_maxg = 1;
    current_w = 1;
}

static void ensure_size(int n) {
    if (n < W_alloc) return;
    int na = W_alloc;
    while (na <= n) na *= 2;
    B_arr = realloc(B_arr, na * sizeof(int));
    g_const_arr = realloc(g_const_arr, na * sizeof(int));
    const_arr = realloc(const_arr, na * sizeof(int));
    G_arr = realloc(G_arr, na * sizeof(int *));
    memset(B_arr + W_alloc, 0, (na - W_alloc) * sizeof(int));
    memset(g_const_arr + W_alloc, 0, (na - W_alloc) * sizeof(int));
    memset(const_arr + W_alloc, 0, (na - W_alloc) * sizeof(int));
    memset(G_arr + W_alloc, 0, (na - W_alloc) * sizeof(int *));
    W_alloc = na;
}

static int getG(int wi, int hi) {
    if (wi > hi) { int t = wi; wi = hi; hi = t; }
    if (wi <= 1) return 0;
    if (hi <= B_arr[wi]) return G_arr[wi][hi];
    return const_arr[wi];
}

/* Mark value v as visited in current stamp, expanding visited if needed. */
static inline void mark_visited(int v) {
    if (v < vis_size) {
        visited[v] = stamp;
    } else {
        int ns = vis_size;
        while (ns <= v) ns *= 2;
        visited = realloc(visited, ns * sizeof(int));
        memset(visited + vis_size, 0, (ns - vis_size) * sizeof(int));
        vis_size = ns;
        visited[v] = stamp;
    }
}

/* Compute nimbers g_w[h] for h in [h_start..h_end]. */
static void compute_h_range(int h_start, int h_end, int *g_w,
                            int **Vs, int *t_as, int half_w) {
    for (int h = h_start; h <= h_end; h++) {
        stamp++;
        int half_h = h >> 1;
        for (int a_idx = 0; a_idx < half_w; a_idx++) {
            int *V = Vs[a_idx];
            int ta = t_as[a_idx];
            int maxb;
            if (half_h >= ta) {
                visited[0] = stamp; /* 0 is reachable from tail-tail */
                maxb = ta - 1;
                if (maxb > half_h) maxb = half_h;
            } else {
                maxb = half_h;
            }
            for (int b = 1; b <= maxb; b++) {
                int v = V[b] ^ V[h - b];
                mark_visited(v);
            }
        }
        int m = 0;
        while (visited[m] == stamp) m++;
        g_w[h] = m;
    }
}

static void compute_next(void) {
    int w = current_w + 1;
    ensure_size(w);

    int limit = (w > 2 * prev_maxg) ? w : 2 * prev_maxg;

    int *g_w = calloc(limit + 1, sizeof(int));

    /* Fill h < w by symmetry */
    int upto = (w < limit + 1) ? w : limit + 1;
    for (int h = 1; h < upto; h++)
        g_w[h] = getG(h, w);

    int half_w = w >> 1;

    /* Precompute V arrays and t_as for each vertical split */
    int **Vs = malloc(half_w * sizeof(int *));
    int *t_as = malloc(half_w * sizeof(int));
    for (int a = 1; a <= half_w; a++) {
        int bw = w - a;
        int ta = (g_const_arr[a] > g_const_arr[bw])
                     ? g_const_arr[a] : g_const_arr[bw];
        t_as[a - 1] = ta;
        int *V = calloc(limit + 1, sizeof(int));
        for (int b = 1; b <= limit; b++)
            V[b] = getG(a, b) ^ getG(bw, b);
        Vs[a - 1] = V;
    }

    compute_h_range(w, limit, g_w, Vs, t_as, half_w);

    /* Track last change index */
    int last = 1;
    for (int h = 2; h <= limit; h++)
        if (g_w[h] != g_w[h - 1]) last = h;

    /* Extend until stability: constant on [last, 2*last] => constant forever */
    while (1) {
        int need = 2 * last;
        if (need <= limit) {
            int v0 = g_w[last];
            int ok = 1;
            for (int hh = last + 1; hh <= need; hh++) {
                if (g_w[hh] != v0) { ok = 0; break; }
            }
            if (ok) break;
        }

        int new_limit = limit * 2;
        if (need > new_limit) new_limit = need;
        if (w > new_limit) new_limit = w;

        g_w = realloc(g_w, (new_limit + 1) * sizeof(int));
        memset(g_w + limit + 1, 0, (new_limit - limit) * sizeof(int));

        for (int a_idx = 0; a_idx < half_w; a_idx++) {
            int a = a_idx + 1;
            int bw = w - a;
            Vs[a_idx] = realloc(Vs[a_idx], (new_limit + 1) * sizeof(int));
            int *V = Vs[a_idx];
            for (int b = limit + 1; b <= new_limit; b++)
                V[b] = getG(a, b) ^ getG(bw, b);
        }

        compute_h_range(limit + 1, new_limit, g_w, Vs, t_as, half_w);
        limit = new_limit;

        for (int h = last + 1; h <= limit; h++)
            if (g_w[h] != g_w[h - 1]) last = h;
    }

    B_arr[w] = limit;
    g_const_arr[w] = last;
    const_arr[w] = g_w[last];
    G_arr[w] = g_w;
    current_w = w;
    if (last > prev_maxg) prev_maxg = last;

    for (int a_idx = 0; a_idx < half_w; a_idx++) free(Vs[a_idx]);
    free(Vs);
    free(t_as);
}

static void compute_upto(int W) {
    while (current_w < W) compute_next();
}

/* ----------------- Counting ----------------- */

/* Directly count winning moves for a single rectangle w x h. */
static long long count_winning_moves_C(int w, int h) {
    long long cnt = 0;
    for (int x = 1; x < w; x++) {
        for (int y = 1; y < h; y++) {
            int v = getG(x, y) ^ getG(w - x, y)
                  ^ getG(x, h - y) ^ getG(w - x, h - y);
            if (v == 0) cnt++;
        }
    }
    return cnt;
}

/* Compute D(W,H) via quadruple counting over ordered splits. */
static long long compute_D(int W, int H) {
    int max_m = 0;
    for (int w = 1; w <= W; w++)
        if (g_const_arr[w] > max_m) max_m = g_const_arr[w];

    int max_k = (max_m < H - 1) ? max_m : H - 1;

    /* Precompute row[u][k] = G(u,k) for u=1..W, k=1..max_k */
    int **row = malloc((W + 1) * sizeof(int *));
    for (int u = 0; u <= W; u++)
        row[u] = calloc(max_k + 1, sizeof(int));
    for (int u = 1; u <= W; u++)
        for (int k = 1; k <= max_k; k++)
            row[u][k] = getG(u, k);

    /* Stamp-based counting arrays for the fast branch */
    int cnts_sz = 65536;
    int *cnts_stamp = calloc(cnts_sz, sizeof(int));
    int *cnts_count = calloc(cnts_sz, sizeof(int));
    int *cnts_list = malloc(cnts_sz * sizeof(int));
    int cur_stamp = 0;

    long long total = 0;

    for (int i = 1; i < W; i++) {
        int *ri = row[i];
        for (int j = 1; j <= W - i; j++) {
            int *rj = row[j];
            int m = (g_const_arr[i] > g_const_arr[j])
                        ? g_const_arr[i] : g_const_arr[j];
            int const_t = const_arr[i] ^ const_arr[j];

            if (H >= 2 * m) {
                /* Fast branch */
                int L = m - 1;
                long long countA = 0, countB = 0, countC = 0;
                int count_eq = 0;
                long long sum_k_eq = 0;

                if (L > 0) {
                    cur_stamp++;
                    int nvisited = 0;
                    for (int k = 1; k <= L; k++) {
                        int t = ri[k] ^ rj[k];
                        if (t >= cnts_sz) {
                            /* Expand (rare) */
                            int ns = cnts_sz;
                            while (ns <= t) ns *= 2;
                            cnts_stamp = realloc(cnts_stamp, ns * sizeof(int));
                            cnts_count = realloc(cnts_count, ns * sizeof(int));
                            cnts_list = realloc(cnts_list, ns * sizeof(int));
                            memset(cnts_stamp + cnts_sz, 0,
                                   (ns - cnts_sz) * sizeof(int));
                            cnts_sz = ns;
                        }
                        if (cnts_stamp[t] != cur_stamp) {
                            cnts_stamp[t] = cur_stamp;
                            cnts_count[t] = 1;
                            cnts_list[nvisited++] = t;
                        } else {
                            cnts_count[t]++;
                        }
                        if (t == const_t) {
                            count_eq++;
                            sum_k_eq += k;
                        }
                    }
                    for (int idx = 0; idx < nvisited; idx++) {
                        int c = cnts_count[cnts_list[idx]];
                        countA += (long long)c * c;
                    }
                }

                if (count_eq) {
                    countB = 2LL * ((long long)count_eq * (H - m + 1)
                                     - sum_k_eq);
                }

                long long S = (long long)H - 2 * m;
                if (S >= 0)
                    countC = (S + 1) * (S + 2) / 2;

                total += countA + countB + countC;
            } else {
                /* Generic smaller-H branch (used for self-tests). */
                if (H > 200000) {
                    fprintf(stderr,
                            "p933: slow branch with large H=%d\n", H);
                    exit(1);
                }
                int *t_arr = malloc(H * sizeof(int));
                for (int k = 1; k < H; k++)
                    t_arr[k] = (k <= max_k) ? ri[k] ^ rj[k] : const_t;
                for (int k = 1; k < H; k++) {
                    int tk = t_arr[k];
                    for (int l = 1; l <= H - k; l++) {
                        if (tk == t_arr[l]) total++;
                    }
                }
                free(t_arr);
            }
        }
    }

    for (int u = 0; u <= W; u++) free(row[u]);
    free(row);
    free(cnts_stamp);
    free(cnts_count);
    free(cnts_list);

    return total;
}

/* ----------------- Entry point ----------------- */

long long p933_native(void) {
    grundy_init();
    compute_upto(123);

    /* Self-tests from the problem statement */
    if (count_winning_moves_C(5, 3) != 4) {
        fprintf(stderr, "p933 self-test C(5,3) failed\n");
        return -1;
    }
    if (compute_D(12, 123) != 327398) {
        fprintf(stderr, "p933 self-test D(12,123) failed\n");
        return -1;
    }

    return compute_D(123, 1234567);
}
