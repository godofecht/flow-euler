#include <stdlib.h>
#include <math.h>
#include <string.h>

static const double SQRT2 = 1.4142135623730950488016887242096980785696718753769;
static const double EPS = 1e-12;

/* ---- dynamic arrays ---- */

typedef struct { double *data; int len; int cap; } DVec;
typedef struct { int    *data; int len; int cap; } IVec;

static void dvec_push(DVec *v, double x) {
    if (v->len >= v->cap) {
        v->cap = v->cap ? v->cap * 2 : 64;
        v->data = realloc(v->data, v->cap * sizeof(double));
    }
    v->data[v->len++] = x;
}

static void ivec_push(IVec *v, int x) {
    if (v->len >= v->cap) {
        v->cap = v->cap ? v->cap * 2 : 64;
        v->data = realloc(v->data, v->cap * sizeof(int));
    }
    v->data[v->len++] = x;
}

static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

/* sort + deduplicate in place, return new length */
static int dedup_double(double *arr, int len) {
    if (len <= 1) return len;
    int j = 0;
    for (int i = 1; i < len; i++) {
        if (arr[i] != arr[j]) arr[++j] = arr[i];
    }
    return j + 1;
}

/* ---- _generate_ring ---- */

static DVec generate_ring(double max_l) {
    DVec vals = {0};
    dvec_push(&vals, 0.0);
    int max_b = (int)(max_l / SQRT2) + 1;
    for (int b = 0; b <= max_b; b++) {
        double base = b * SQRT2;
        int max_a = (int)(max_l - base + 1e-12);
        for (int a = 0; a <= max_a; a++)
            dvec_push(&vals, (double)a + base);
    }
    qsort(vals.data, vals.len, sizeof(double), cmp_double);
    vals.len = dedup_double(vals.data, vals.len);
    return vals;
}

/* ---- _compute_grundy_intervals ---- */

static int upper_index(const double *arr, int len, double x) {
    int lo = 0, hi = len;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (arr[mid] <= x) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

typedef struct {
    DVec starts, ends;
    IVec grundy;
} GrundyResult;

static GrundyResult compute_grundy_intervals(double max_l) {
    DVec vals = generate_ring(max_l);
    int n = vals.len - 1;

    DVec starts = {0}, ends = {0};
    IVec grundy = {0};

    int moves_cap = 256;
    char *moves = calloc(moves_cap, sizeof(char));
    int max_g = 0;

    for (int idx = 0; idx < n; idx++) {
        double a = vals.data[idx];
        double b = vals.data[idx + 1];
        double L = (a + b) * 0.5;

        for (int i = 0; i <= max_g; i++) moves[i] = 0;

        double xs[2] = {1.0, SQRT2};
        for (int xi = 0; xi < 2; xi++) {
            double x = xs[xi];
            if (L < x || starts.len == 0) continue;
            double S = L - x;
            int u_idx = upper_index(starts.data, starts.len, S) - 1;
            int t_idx = 0;
            while (t_idx <= u_idx) {
                double t_start = starts.data[t_idx];
                if (t_start >= S) break;
                double t_end = ends.data[t_idx];
                if (t_end > S) t_end = S;
                double u_start = starts.data[u_idx];
                double u_end = ends.data[u_idx];
                if (u_end > S) u_end = S;
                double left = S - u_end;
                if (t_start > left) left = t_start;
                double right = S - u_start;
                if (t_end < right) right = t_end;
                if (left < right) {
                    int g = grundy.data[t_idx] ^ grundy.data[u_idx];
                    if (g >= moves_cap) {
                        int new_cap = moves_cap;
                        while (g >= new_cap) new_cap *= 2;
                        moves = realloc(moves, new_cap * sizeof(char));
                        memset(moves + moves_cap, 0,
                               (new_cap - moves_cap) * sizeof(char));
                        moves_cap = new_cap;
                    }
                    moves[g] = 1;
                    if (g > max_g) max_g = g;
                }
                if (t_end < S - u_start) t_idx++;
                else u_idx--;
            }
        }

        int g = 0;
        while (g < moves_cap && moves[g]) g++;
        if (g >= moves_cap) {
            int new_cap = moves_cap * 2;
            moves = realloc(moves, new_cap * sizeof(char));
            memset(moves + moves_cap, 0,
                   (new_cap - moves_cap) * sizeof(char));
            moves_cap = new_cap;
        }

        if (grundy.len > 0 &&
            grundy.data[grundy.len - 1] == g &&
            fabs(ends.data[ends.len - 1] - a) < EPS) {
            ends.data[ends.len - 1] = b;
        } else {
            dvec_push(&starts, a);
            dvec_push(&ends, b);
            ivec_push(&grundy, g);
        }
    }

    free(moves);
    free(vals.data);

    GrundyResult r = { starts, ends, grundy };
    return r;
}

/* ---- _build_w_segments ---- */

typedef struct { double pos; int delta; } Event;

static int cmp_event(const void *a, const void *b) {
    const Event *ea = a, *eb = b;
    if (ea->pos < eb->pos) return -1;
    if (ea->pos > eb->pos) return 1;
    return 0;
}

typedef struct {
    DVec seg_starts, seg_ends, seg_slopes, seg_vals;
} SegData;

static SegData build_w_segments(DVec starts, DVec ends, IVec grundy,
                                double max_s) {
    int max_g = 0;
    for (int i = 0; i < grundy.len; i++)
        if (grundy.data[i] > max_g) max_g = grundy.data[i];

    /* group intervals by grundy number */
    int *gsz = calloc(max_g + 1, sizeof(int));
    for (int i = 0; i < grundy.len; i++) gsz[grundy.data[i]]++;

    double **ga = calloc(max_g + 1, sizeof(double *));
    double **gb = calloc(max_g + 1, sizeof(double *));
    int   *gidx = calloc(max_g + 1, sizeof(int));
    for (int g = 0; g <= max_g; g++) {
        ga[g] = malloc(gsz[g] * sizeof(double));
        gb[g] = malloc(gsz[g] * sizeof(double));
    }
    for (int i = 0; i < grundy.len; i++) {
        int g = grundy.data[i];
        ga[g][gidx[g]] = starts.data[i];
        gb[g][gidx[g]] = ends.data[i];
        gidx[g]++;
    }

    /* generate events */
    Event *events = NULL;
    int ev_count = 0, ev_cap = 0;

    for (int g = 0; g <= max_g; g++) {
        int m = gsz[g];
        for (int i = 0; i < m; i++) {
            double a1 = ga[g][i], b1 = gb[g][i];
            for (int j = i; j < m; j++) {
                double a2 = ga[g][j], b2 = gb[g][j];
                int w = (i == j) ? 1 : 2;
                double p0 = a1 + a2, p1 = a1 + b2, p2 = b1 + a2, p3 = b1 + b2;
                if (p0 > max_s + EPS) break;
                if (p3 < 0.0) continue;
                if (p0 < 0.0) p0 = 0.0;
                if (p3 > max_s) p3 = max_s;
                double q1, q2;
                if (p1 < p2) { q1 = p1; q2 = p2; } else { q1 = p2; q2 = p1; }

                double posv[4]  = { p0, q1, q2, p3 };
                int    deltav[4] = { w,  -w,  -w,  w  };
                for (int k = 0; k < 4; k++) {
                    if (ev_count >= ev_cap) {
                        ev_cap = ev_cap ? ev_cap * 2 : 4096;
                        events = realloc(events, ev_cap * sizeof(Event));
                    }
                    events[ev_count].pos   = posv[k];
                    events[ev_count].delta = deltav[k];
                    ev_count++;
                }
            }
        }
    }

    qsort(events, ev_count, sizeof(Event), cmp_event);

    /* merge events by position (within EPS) */
    Event *merged = NULL;
    int m_count = 0, m_cap = 0;
    if (ev_count > 0) {
        double cur_pos = events[0].pos;
        int    cur_delta = events[0].delta;
        for (int i = 1; i < ev_count; i++) {
            if (fabs(events[i].pos - cur_pos) < EPS) {
                cur_delta += events[i].delta;
            } else {
                if (m_count >= m_cap) {
                    m_cap = m_cap ? m_cap * 2 : 4096;
                    merged = realloc(merged, m_cap * sizeof(Event));
                }
                merged[m_count].pos   = cur_pos;
                merged[m_count].delta = cur_delta;
                m_count++;
                cur_pos   = events[i].pos;
                cur_delta = events[i].delta;
            }
        }
        if (m_count >= m_cap) {
            m_cap = m_cap ? m_cap * 2 : 4096;
            merged = realloc(merged, m_cap * sizeof(Event));
        }
        merged[m_count].pos   = cur_pos;
        merged[m_count].delta = cur_delta;
        m_count++;
    }

    /* build piecewise-linear segments */
    SegData sd = {0};
    double slope = 0.0, val = 0.0, prev = 0.0;
    for (int i = 0; i < m_count; i++) {
        double pos = merged[i].pos;
        int    delta = merged[i].delta;
        if (pos > max_s) break;
        if (pos > prev) {
            dvec_push(&sd.seg_starts, prev);
            dvec_push(&sd.seg_ends,   pos);
            dvec_push(&sd.seg_slopes, slope);
            dvec_push(&sd.seg_vals,   val);
            val += slope * (pos - prev);
            prev = pos;
        }
        slope += delta;
    }
    if (prev < max_s) {
        dvec_push(&sd.seg_starts, prev);
        dvec_push(&sd.seg_ends,   max_s);
        dvec_push(&sd.seg_slopes, slope);
        dvec_push(&sd.seg_vals,   val);
    }

    for (int g = 0; g <= max_g; g++) { free(ga[g]); free(gb[g]); }
    free(ga); free(gb); free(gsz); free(gidx);
    free(events); free(merged);

    return sd;
}

/* ---- _w_value ---- */

static void w_value(const SegData *sd, double x,
                    double *out_val, double *out_slope) {
    int lo = 0, hi = sd->seg_starts.len;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (sd->seg_ends.data[mid] <= x) lo = mid + 1;
        else hi = mid;
    }
    if (lo >= sd->seg_starts.len) {
        *out_val   = 0.0;
        *out_slope = 0.0;
        return;
    }
    double start = sd->seg_starts.data[lo];
    double slope = sd->seg_slopes.data[lo];
    double val   = sd->seg_vals.data[lo] + slope * (x - start);
    *out_val   = val;
    *out_slope = slope;
}

/* ---- _e_value ---- */

static double e_value(double L, const SegData *sd) {
    double w1, s1, w2, s2;
    w_value(sd, L - 1.0,  &w1, &s1);
    w_value(sd, L - SQRT2, &w2, &s2);
    return 0.5 * L * (w1 / (L - 1.0) + w2 / (L - SQRT2));
}

/* ---- _f_value helpers (nested closures lifted to top level) ---- */

static double term_deriv(double m, double b0, double c, double L) {
    return (m * L * L - 2.0 * m * c * L - b0 * c) / ((L - c) * (L - c));
}

static double e_local_func(double L,
                           double m1, double b1,
                           double m2, double b2) {
    return 0.5 * L * ((m1 * L + b1) / (L - 1.0) +
                      (m2 * L + b2) / (L - SQRT2));
}

static double de_local_func(double L,
                            double m1, double b1,
                            double m2, double b2) {
    return 0.5 * (term_deriv(m1, b1, 1.0, L) +
                  term_deriv(m2, b2, SQRT2, L));
}

static double bisect_func(double lo, double hi, double dlo,
                          double m1, double b1,
                          double m2, double b2) {
    for (int iter = 0; iter < 60; iter++) {
        double m = (lo + hi) * 0.5;
        double dm = de_local_func(m, m1, b1, m2, b2);
        if (dm == 0.0) return m;
        if (dm * dlo > 0.0) { lo = m; dlo = dm; }
        else hi = m;
    }
    return (lo + hi) * 0.5;
}

/* ---- _f_value ---- */

static double f_value(double a, double b, const SegData *sd) {
    DVec points = {0};
    dvec_push(&points, a);
    dvec_push(&points, b);
    for (int i = 0; i < sd->seg_starts.len; i++) {
        double p = sd->seg_starts.data[i];
        double v = p + 1.0;
        if (a < v && v < b) dvec_push(&points, v);
        v = p + SQRT2;
        if (a < v && v < b) dvec_push(&points, v);
    }
    for (int i = 0; i < sd->seg_ends.len; i++) {
        double p = sd->seg_ends.data[i];
        double v = p + 1.0;
        if (a < v && v < b) dvec_push(&points, v);
        v = p + SQRT2;
        if (a < v && v < b) dvec_push(&points, v);
    }

    qsort(points.data, points.len, sizeof(double), cmp_double);
    points.len = dedup_double(points.data, points.len);

    double best = -1.0;

    for (int i = 0; i < points.len - 1; i++) {
        double L0 = points.data[i];
        double L1 = points.data[i + 1];
        if (L1 - L0 < 1e-12) continue;
        double mid = (L0 + L1) * 0.5;

        double w1, m1, w2, m2;
        w_value(sd, mid - 1.0,  &w1, &m1);
        w_value(sd, mid - SQRT2, &w2, &m2);
        double b1 = w1 - m1 * mid;
        double b2 = w2 - m2 * mid;

        double checkL[3] = { L0, mid, L1 };
        for (int k = 0; k < 3; k++) {
            double val = e_local_func(checkL[k], m1, b1, m2, b2);
            if (val > best) best = val;
        }

        double left  = L0 + 1e-10;
        double right = L1 - 1e-10;
        if (left >= right) continue;
        double dl = de_local_func(left, m1, b1, m2, b2);
        double dm = de_local_func(mid,  m1, b1, m2, b2);
        double dr = de_local_func(right, m1, b1, m2, b2);

        if (dl * dm < 0.0) {
            double root = bisect_func(left, mid, dl, m1, b1, m2, b2);
            double val  = e_local_func(root, m1, b1, m2, b2);
            if (val > best) best = val;
        }
        if (dm * dr < 0.0) {
            double root = bisect_func(mid, right, dm, m1, b1, m2, b2);
            double val  = e_local_func(root, m1, b1, m2, b2);
            if (val > best) best = val;
        }
    }

    free(points.data);
    return best;
}

/* ---- entry point ---- */

double p644_native(void) {
    double max_l = 500.0;
    GrundyResult gr = compute_grundy_intervals(max_l);
    SegData sd = build_w_segments(gr.starts, gr.ends, gr.grundy, max_l);

    double answer = f_value(200.0, 500.0, &sd);

    free(gr.starts.data);  free(gr.ends.data);  free(gr.grundy.data);
    free(sd.seg_starts.data); free(sd.seg_ends.data);
    free(sd.seg_slopes.data); free(sd.seg_vals.data);

    return answer;
}
