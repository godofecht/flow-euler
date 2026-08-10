
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static int cmp_double(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static int bisect_left(const double *a, int n, double x) {
    int lo = 0, hi = n;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (a[mid] < x) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

static int findp(int *parent, int i) {
    while (parent[i] != i) {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    return i;
}

typedef struct { double *ends; double *vals; int n; } Seg;

static Seg build_piecewise(double l, double g) {
    double domain_end = 1.0 - g;
    int K = (int)(1.5 / g); if (K < 10) K = 10;
    for (;;) {
        if (K > 4000000) { fprintf(stderr, "K overflow\n"); abort(); }
        double *starts = malloc((size_t)K * sizeof(double));
        double *ends = malloc((size_t)K * sizeof(double));
        int *labels = malloc((size_t)K * sizeof(int));
        int nint = 0;
        double *endpoints = malloc((size_t)(2 * K + 2) * sizeof(double));
        int nep = 0;
        endpoints[nep++] = 0.0;
        endpoints[nep++] = domain_end;
        double x = 0.0;
        for (int k = 1; k <= K; k++) {
            x += l;
            x -= (double)(long long)x; /* toward zero, x>=0 */
            double s = x - g, e = x;
            if (e <= 0.0 || s >= domain_end) continue;
            if (s < 0.0) s = 0.0;
            if (e > domain_end) e = domain_end;
            if (s < e) {
                starts[nint] = s; ends[nint] = e; labels[nint] = k;
                endpoints[nep++] = s; endpoints[nep++] = e;
                nint++;
            }
        }
        qsort(endpoints, (size_t)nep, sizeof(double), cmp_double);
        int nu = 0;
        for (int i = 0; i < nep; i++) {
            if (nu == 0 || endpoints[i] != endpoints[nu - 1])
                endpoints[nu++] = endpoints[i];
        }
        int num_cells = nu - 1;
        int *parent = malloc((size_t)(num_cells + 1) * sizeof(int));
        unsigned *values = calloc((size_t)num_cells, sizeof(unsigned));
        for (int i = 0; i <= num_cells; i++) parent[i] = i;
        for (int t = 0; t < nint; t++) {
            int i = bisect_left(endpoints, nu, starts[t]);
            int j = bisect_left(endpoints, nu, ends[t]);
            int idx = findp(parent, i);
            while (idx < j) {
                values[idx] = (unsigned)labels[t];
                parent[idx] = idx + 1;
                idx = findp(parent, idx);
            }
        }
        if (findp(parent, 0) == num_cells) {
            Seg seg = {0};
            seg.ends = malloc((size_t)(num_cells + 1) * sizeof(double));
            seg.vals = malloc((size_t)(num_cells + 1) * sizeof(double));
            unsigned curr = values[0];
            for (int i = 1; i < num_cells; i++) {
                if (values[i] != curr) {
                    seg.ends[seg.n] = endpoints[i];
                    seg.vals[seg.n] = (double)curr * l;
                    seg.n++;
                    curr = values[i];
                }
            }
            seg.ends[seg.n] = domain_end;
            seg.vals[seg.n] = (double)curr * l;
            seg.n++;
            free(starts); free(ends); free(labels); free(endpoints); free(parent); free(values);
            return seg;
        }
        free(starts); free(ends); free(labels); free(endpoints); free(parent); free(values);
        K *= 2;
    }
}

typedef struct { double key; int id; } HI;
static void hpush(HI *h, int *n, HI v) {
    int i = (*n)++;
    h[i] = v;
    while (i) {
        int p = (i - 1) / 2;
        if (h[p].key <= h[i].key) break;
        HI t = h[p]; h[p] = h[i]; h[i] = t; i = p;
    }
}
static HI hpop(HI *h, int *n) {
    HI r = h[0];
    h[0] = h[--(*n)];
    int i = 0;
    for (;;) {
        int l = 2 * i + 1, rg = 2 * i + 2, sm = i;
        if (l < *n && h[l].key < h[sm].key) sm = l;
        if (rg < *n && h[rg].key < h[sm].key) sm = rg;
        if (sm == i) break;
        HI t = h[i]; h[i] = h[sm]; h[sm] = t; i = sm;
    }
    return r;
}

static double merge_max_sum(Seg *segs, int P, double domain_end) {
    int *idx = calloc((size_t)P, sizeof(int));
    double *cur = malloc((size_t)P * sizeof(double));
    double total = 0.0;
    for (int i = 0; i < P; i++) { cur[i] = segs[i].vals[0]; total += cur[i]; }
    HI *heap = malloc((size_t)(P * 2 + 8) * sizeof(HI));
    int hn = 0;
    for (int i = 0; i < P; i++) { HI it = {segs[i].ends[0], i}; hpush(heap, &hn, it); }
    double cur_pos = 0.0, best = total;
    const double eps = 1e-15;
    while (hn > 0) {
        double boundary = heap[0].key;
        if (boundary > cur_pos + 1e-18 && total > best) best = total;
        int aff[512]; int na = 0;
        while (hn > 0 && fabs(heap[0].key - boundary) <= eps) {
            if (na < 512) aff[na++] = hpop(heap, &hn).id;
            else hpop(heap, &hn);
        }
        cur_pos = boundary;
        if (cur_pos >= domain_end - 1e-15) break;
        for (int a = 0; a < na; a++) {
            int i = aff[a];
            double old = cur[i];
            idx[i]++;
            double nv = segs[i].vals[idx[i]];
            cur[i] = nv;
            total += nv - old;
            HI it = {segs[i].ends[idx[i]], i};
            hpush(heap, &hn, it);
        }
    }
    free(idx); free(cur); free(heap);
    return best;
}

static int *primes_up_to(int n, int *outc) {
    char *sv = calloc((size_t)n + 1, 1);
    memset(sv, 1, (size_t)n + 1); sv[0] = sv[1] = 0;
    int r = (int)sqrt((double)n);
    for (int i = 2; i <= r; i++) if (sv[i])
        for (long j = (long)i * i; j <= n; j += i) sv[j] = 0;
    int c = 0; for (int i = 2; i <= n; i++) if (sv[i]) c++;
    int *ps = malloc((size_t)c * sizeof(int));
    c = 0; for (int i = 2; i <= n; i++) if (sv[i]) ps[c++] = i;
    free(sv); *outc = c; return ps;
}

static double M(int n, double g) {
    int npc; int *ps = primes_up_to(n, &npc);
    Seg *segs = malloc((size_t)npc * sizeof(Seg));
    for (int i = 0; i < npc; i++) segs[i] = build_piecewise(1.0 / sqrt((double)ps[i]), g);
    double ans = merge_max_sum(segs, npc, 1.0 - g);
    for (int i = 0; i < npc; i++) { free(segs[i].ends); free(segs[i].vals); }
    free(segs); free(ps);
    return ans;
}

double pe576_answer(void) { return M(100, 0.00002); }
