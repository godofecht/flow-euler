/* Project Euler 742: Minimum Area of a Convex Grid Polygon
   Port of the Python reference solver. */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct { int a, b; } Pair;

static int gcd(int a, int b) {
    while (b) { int t = a % b; a = b; b = t; }
    return a;
}

/* primitive pairs (a,b) with 1<=a,b<=limit, gcd(a,b)==1 */
static Pair *g_pairs = NULL;
static int g_npairs = 0;
static int g_pair_cap = 0;

static void primitive_pairs(int limit) {
    if (g_pair_cap < limit * limit) {
        free(g_pairs);
        g_pair_cap = limit * limit;
        g_pairs = (Pair*)malloc(g_pair_cap * sizeof(Pair));
    }
    g_npairs = 0;
    for (int a = 1; a <= limit; a++) {
        for (int b = 1; b <= limit; b++) {
            if (gcd(a, b) == 1) {
                g_pairs[g_npairs].a = a;
                g_pairs[g_npairs].b = b;
                g_npairs++;
            }
        }
    }
}

/* Select k pairs minimizing a^2 + t*b^2, with tie-break on a+b then a then b.
   Uses a max-heap of size k. */
static void n_smallest(Pair *out, int k, double t) {
    /* heap entries: weight, tie1(a+b), a, b. Max-heap by (weight, tie1, a, b). */
    typedef struct { double w; int tie1, a, b; } HEntry;
    HEntry *heap = (HEntry*)malloc(k * sizeof(HEntry));
    int hs = 0;

    for (int idx = 0; idx < g_npairs; idx++) {
        int a = g_pairs[idx].a, b = g_pairs[idx].b;
        double w = (double)(a*a) + t * (double)(b*b);
        int tie1 = a + b;
        if (hs < k) {
            /* push */
            int i = hs++;
            heap[i].w = w; heap[i].tie1 = tie1; heap[i].a = a; heap[i].b = b;
            /* sift up */
            while (i > 0) {
                int p = (i - 1) / 2;
                /* max-heap: parent should be >= child */
                if (heap[p].w < heap[i].w ||
                    (heap[p].w == heap[i].w && heap[p].tie1 < heap[i].tie1) ||
                    (heap[p].w == heap[i].w && heap[p].tie1 == heap[i].tie1 && heap[p].a < heap[i].a) ||
                    (heap[p].w == heap[i].w && heap[p].tie1 == heap[i].tie1 && heap[p].a == heap[i].a && heap[p].b < heap[i].b)) {
                    HEntry tmp = heap[p]; heap[p] = heap[i]; heap[i] = tmp;
                    i = p;
                } else break;
            }
        } else {
            /* compare with root (max) */
            int better = (w < heap[0].w) ||
                (w == heap[0].w && tie1 < heap[0].tie1) ||
                (w == heap[0].w && tie1 == heap[0].tie1 && a < heap[0].a) ||
                (w == heap[0].w && tie1 == heap[0].tie1 && a == heap[0].a && b < heap[0].b);
            if (better) {
                /* replace root */
                heap[0].w = w; heap[0].tie1 = tie1; heap[0].a = a; heap[0].b = b;
                /* sift down */
                int i = 0;
                while (1) {
                    int l = 2*i+1, r = 2*i+2, largest = i;
                    if (l < hs && (heap[l].w > heap[largest].w ||
                        (heap[l].w == heap[largest].w && heap[l].tie1 > heap[largest].tie1) ||
                        (heap[l].w == heap[largest].w && heap[l].tie1 == heap[largest].tie1 && heap[l].a > heap[largest].a) ||
                        (heap[l].w == heap[largest].w && heap[l].tie1 == heap[largest].tie1 && heap[l].a == heap[largest].a && heap[l].b > heap[largest].b)))
                        largest = l;
                    if (r < hs && (heap[r].w > heap[largest].w ||
                        (heap[r].w == heap[largest].w && heap[r].tie1 > heap[largest].tie1) ||
                        (heap[r].w == heap[largest].w && heap[r].tie1 == heap[largest].tie1 && heap[r].a > heap[largest].a) ||
                        (heap[r].w == heap[largest].w && heap[r].tie1 == heap[largest].tie1 && heap[r].a == heap[largest].a && heap[r].b > heap[largest].b)))
                        largest = r;
                    if (largest != i) {
                        HEntry tmp = heap[i]; heap[i] = heap[largest]; heap[largest] = tmp;
                        i = largest;
                    } else break;
                }
            }
        }
    }

    for (int i = 0; i < k; i++) {
        out[i].a = heap[i].a;
        out[i].b = heap[i].b;
    }
    free(heap);
}

static int cmp_slope(const void *pa, const void *pb) {
    const Pair *a = (const Pair*)pa;
    const Pair *b = (const Pair*)pb;
    /* sort by slope b/a ascending, then a, then b */
    double sa = (double)a->b / (double)a->a;
    double sb = (double)b->b / (double)b->a;
    if (sa < sb) return -1;
    if (sa > sb) return 1;
    if (a->a < b->a) return -1;
    if (a->a > b->a) return 1;
    if (a->b < b->b) return -1;
    if (a->b > b->b) return 1;
    return 0;
}

static long long area_from_half_edges(Pair *half, int m) {
    long long px = 0, py = 0, area = 0;
    for (int i = 0; i < m; i++) {
        area += px * half[i].b - py * half[i].a;
        px += half[i].a;
        py += half[i].b;
    }
    return area;
}

static long long polygon_area_from_interior(Pair *interior, int k) {
    /* sort interior by slope */
    qsort(interior, k, sizeof(Pair), cmp_slope);

    /* half = (1,0) + interior_sorted + (0,1) + mirrored reversed interior */
    int m = 1 + k + 1 + k;
    Pair *half = (Pair*)malloc(m * sizeof(Pair));
    int idx = 0;
    half[idx].a = 1; half[idx].b = 0; idx++;
    for (int i = 0; i < k; i++) { half[idx] = interior[i]; idx++; }
    half[idx].a = 0; half[idx].b = 1; idx++;
    for (int i = k - 1; i >= 0; i--) { half[idx].a = -interior[i].a; half[idx].b = interior[i].b; idx++; }

    long long area = area_from_half_edges(half, m);
    free(half);
    return area;
}

static long long compute_A(int N) {
    if (N < 4 || N % 4 != 0) return 0;
    int k = (N - 4) / 4;
    if (k == 0) return 1;

    int limit = 40;
    primitive_pairs(limit);

    long long best_area = -1;
    Pair *chosen = (Pair*)malloc(k * sizeof(Pair));

    for (int tn = 1; tn <= 1000; tn++) {
        double t = (double)tn / 1000.0;

        /* ensure enough candidates */
        while (1) {
            n_smallest(chosen, k, t);
            int max_a = 0, max_b = 0;
            for (int i = 0; i < k; i++) {
                if (chosen[i].a > max_a) max_a = chosen[i].a;
                if (chosen[i].b > max_b) max_b = chosen[i].b;
            }
            if (max_a < limit && max_b < limit) break;
            limit *= 2;
            primitive_pairs(limit);
        }

        long long area = polygon_area_from_interior(chosen, k);
        if (best_area == -1 || area < best_area) best_area = area;
    }

    free(chosen);
    return best_area;
}

long long p742_native(void) {
    return compute_A(1000);
}
