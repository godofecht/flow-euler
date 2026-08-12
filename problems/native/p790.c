/*
 * Project Euler 790 - Clock Grid
 *
 * Compute C(10^5) for a 50515093 x 50515093 grid where each of 10^5 rectangles
 * increments the "hour" (overlap count mod 12) of covered cells.
 *
 * Algorithm: sweep line over x with a segment tree on compressed y-coordinates.
 * Each segment tree node tracks 12 buckets (total y-length per overlap residue
 * mod 12).  Range updates rotate the buckets; the sweep accumulates grid-point
 * counts per residue, then converts residues to displayed hours.
 *
 * Ported from the Python reference solver.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t  i64;
typedef uint64_t u64;
typedef uint8_t  u8;

#define M_VAL 50515093LL
#define S0    290797LL

/* ---- comparators ---- */
static int cmp_i64(const void *a, const void *b) {
    i64 va = *(const i64 *)a, vb = *(const i64 *)b;
    return va < vb ? -1 : va > vb ? 1 : 0;
}

typedef struct {
    i64  x;
    int  shift;
    i64  yl;    /* y value, later overwritten with y-index */
    i64  yh1;   /* y value, later overwritten with y-index */
} Event;

static int cmp_event_x(const void *a, const void *b) {
    const Event *ea = a, *eb = b;
    return ea->x < eb->x ? -1 : ea->x > eb->x ? 1 : 0;
}

/* ---- binary search: index of val in sorted arr (val guaranteed present) ---- */
static int bsearch_idx(const i64 *arr, int n, i64 val) {
    int lo = 0, hi = n - 1;
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (arr[mid] < val) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

/* ---- segment tree (12-bucket lazy rotation) ---- */
static u64 *seg;
static u8  *lazy;
static i64 *g_y_vals;

static void build_tree(int node, int l, int r) {
    int base = node * 12;
    if (r - l == 1) {
        seg[base] = (u64)(g_y_vals[l + 1] - g_y_vals[l]);
        return;
    }
    int mid = (l + r) >> 1;
    int left = node << 1;
    build_tree(left, l, mid);
    build_tree(left + 1, mid, r);
    seg[base] = seg[left * 12] + seg[(left + 1) * 12];
}

static void apply_rot(int node, int shift) {
    if (shift == 0) return;
    shift %= 12;
    int base = node * 12;
    u64 tmp[12];
    for (int i = 0; i < 12; i++)
        tmp[i] = seg[base + ((i - shift + 12) % 12)];
    for (int i = 0; i < 12; i++)
        seg[base + i] = tmp[i];
    int v = lazy[node] + shift;
    lazy[node] = (u8)(v >= 12 ? v - 12 : v);
}

static void push_down(int node) {
    int s = lazy[node];
    if (s) {
        int left = node << 1;
        apply_rot(left, s);
        apply_rot(left + 1, s);
        lazy[node] = 0;
    }
}

static void pull_up(int node) {
    int base = node * 12;
    int bl = (node << 1) * 12;
    int br = ((node << 1) + 1) * 12;
    for (int i = 0; i < 12; i++)
        seg[base + i] = seg[bl + i] + seg[br + i];
}

static void update_tree(int node, int l, int r, int ql, int qr, int shift) {
    if (ql <= l && r <= qr) {
        apply_rot(node, shift);
        return;
    }
    push_down(node);
    int mid = (l + r) >> 1;
    int left = node << 1;
    if (ql < mid) update_tree(left, l, mid, ql, qr, shift);
    if (qr > mid) update_tree(left + 1, mid, r, ql, qr, shift);
    pull_up(node);
}

long long p790_native(void) {
    int t = 100000;

    /* ---- generate rectangles, collect x/y coords, build events ---- */
    int cap = 2 + 2 * t;
    i64 *x_arr = malloc((size_t)cap * sizeof(i64));
    i64 *y_arr = malloc((size_t)cap * sizeof(i64));
    int nx = 0, ny = 0;

    x_arr[nx++] = 0;  x_arr[nx++] = M_VAL;
    y_arr[ny++] = 0;  y_arr[ny++] = M_VAL;

    Event *events = malloc((size_t)2 * t * sizeof(Event));
    int nev = 0;

    i64 s = S0;
    for (int i = 0; i < t; i++) {
        i64 x1 = s;  s = (s * s) % M_VAL;
        i64 x2 = s;  s = (s * s) % M_VAL;
        i64 y1 = s;  s = (s * s) % M_VAL;
        i64 y2 = s;  s = (s * s) % M_VAL;

        i64 xl, xh, yl, yh;
        if (x1 <= x2) { xl = x1; xh = x2; }
        else          { xl = x2; xh = x1; }
        if (y1 <= y2) { yl = y1; yh = y2; }
        else          { yl = y2; yh = y1; }

        i64 xh1 = xh + 1;
        i64 yh1 = yh + 1;

        x_arr[nx++] = xl;  x_arr[nx++] = xh1;
        y_arr[ny++] = yl;  y_arr[ny++] = yh1;

        events[nev].x = xl;   events[nev].shift = 1;
        events[nev].yl = yl;  events[nev].yh1 = yh1;
        nev++;

        events[nev].x = xh1;   events[nev].shift = 11;
        events[nev].yl = yl;   events[nev].yh1 = yh1;
        nev++;
    }

    /* ---- sort and deduplicate x and y ---- */
    qsort(x_arr, nx, sizeof(i64), cmp_i64);
    qsort(y_arr, ny, sizeof(i64), cmp_i64);

    int num_x = 1;
    for (int i = 1; i < nx; i++)
        if (x_arr[i] != x_arr[num_x - 1])
            x_arr[num_x++] = x_arr[i];

    int num_y = 1;
    for (int i = 1; i < ny; i++)
        if (y_arr[i] != y_arr[num_y - 1])
            y_arr[num_y++] = y_arr[i];

    /* ---- map event y-values to compressed indices ---- */
    for (int i = 0; i < nev; i++) {
        events[i].yl  = bsearch_idx(y_arr, num_y, events[i].yl);
        events[i].yh1 = bsearch_idx(y_arr, num_y, events[i].yh1);
    }

    /* ---- sort events by x ---- */
    qsort(events, nev, sizeof(Event), cmp_event_x);

    /* ---- build segment tree over y-intervals ---- */
    int m = num_y - 1;
    int tree_size = 4 * m + 5;
    seg  = calloc((size_t)12 * tree_size, sizeof(u64));
    lazy = calloc((size_t)tree_size, 1);
    g_y_vals = y_arr;

    build_tree(1, 0, m);

    /* ---- sweep along x ---- */
    i64 counts[12];
    memset(counts, 0, sizeof(counts));

    int ev_i = 0;
    for (int i = 0; i < num_x - 1; i++) {
        i64 x = x_arr[i];
        while (ev_i < nev && events[ev_i].x == x) {
            update_tree(1, 0, m,
                        (int)events[ev_i].yl, (int)events[ev_i].yh1,
                        events[ev_i].shift);
            ev_i++;
        }

        i64 width = x_arr[i + 1] - x;
        if (width) {
            int base = 12;  /* node 1 */
            for (int r = 0; r < 12; r++)
                counts[r] += width * (i64)seg[base + r];
        }
    }

    /* ---- convert overlap residues to displayed hours ---- */
    i64 total = 12 * counts[0];
    for (int r = 1; r < 12; r++)
        total += (i64)r * counts[r];

    free(seg);  free(lazy);
    free(x_arr); free(y_arr);
    free(events);

    return total;
}
