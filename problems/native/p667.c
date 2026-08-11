/* Project Euler 667: Moving Pentagon
 *
 * Port of the Python reference solver.
 * Nested golden-section + bisection optimization over pentagon shape and scale.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct { double x, y; } Point;

static double heron(double a, double b, double c) {
    double s = (a + b + c) * 0.5;
    double v = s * (s - a) * (s - b) * (s - c);
    return v > 0.0 ? sqrt(v) : 0.0;
}

static void rotate_points(const Point *in, Point *out, int n, double theta) {
    double c = cos(theta), s = sin(theta);
    for (int i = 0; i < n; i++) {
        out[i].x = c * in[i].x - s * in[i].y;
        out[i].y = s * in[i].x + c * in[i].y;
    }
}

static double min_x_above_y(const Point *pts, int n, double ythr) {
    double mn = INFINITY;
    for (int i = 0; i < n; i++) {
        const Point *p1 = &pts[i];
        const Point *p2 = &pts[(i + 1) % n];
        if (p1->y >= ythr && p1->x < mn) mn = p1->x;
        if (p2->y >= ythr && p2->x < mn) mn = p2->x;
        double dy = p2->y - p1->y;
        if (dy == 0.0) continue;
        double t = (ythr - p1->y) / dy;
        if (t >= 0.0 && t <= 1.0) {
            double x = p1->x + (p2->x - p1->x) * t;
            if (x < mn) mn = x;
        }
    }
    return mn;
}

static int build_unit_pentagon(double r, Point *out) {
    if (!(r > 0.5 && r < 2.0)) return 0;
    double h2 = r * r - 0.25;
    if (h2 <= 0.0) return 0;
    double h = sqrt(h2);
    Point A = {0.0, 0.0};
    Point E = {1.0, 0.0};
    Point C = {0.5, h};
    double d = r;
    if (d >= 2.0) return 0;
    double k2 = 1.0 - (d * 0.5) * (d * 0.5);
    if (k2 <= 0.0) return 0;
    double k = sqrt(k2);
    double Mx = (A.x + C.x) * 0.5;
    double My = (A.y + C.y) * 0.5;
    double ux = -(C.y - A.y) / d;
    double uy = (C.x - A.x) / d;
    Point B = {Mx + k * ux, My + k * uy};
    Point D = {1.0 - B.x, B.y};
    out[0] = A; out[1] = B; out[2] = C; out[3] = D; out[4] = E;
    return 1;
}

static double base_area(double r) {
    return 2.0 * heron(1.0, 1.0, r) + heron(r, r, 1.0);
}

static double clearance_for_theta(const Point *pts_rot, double min_y, double max_x,
                                   double scale, double eps_y, int n) {
    double ythr = min_y + (1.0 + eps_y) / scale;
    double x_min = min_x_above_y(pts_rot, n, ythr);
    if (isinf(x_min)) return INFINITY;
    return 1.0 + scale * (x_min - max_x);
}

typedef struct {
    Point pts[5];
    double miny, maxx;
} RotData;

static double min_clearance(const Point *points, RotData *precomp, double *thetas,
                            int n_theta, double scale, double eps_y, int local_k, int local_iters) {
    double *vals = malloc(n_theta * sizeof(double));
    double best = INFINITY;
    for (int i = 0; i < n_theta; i++) {
        vals[i] = clearance_for_theta(precomp[i].pts, precomp[i].miny, precomp[i].maxx, scale, eps_y, 5);
        if (vals[i] < best) best = vals[i];
    }

    /* Sort indices by vals */
    int *idx_sorted = malloc(n_theta * sizeof(int));
    for (int i = 0; i < n_theta; i++) idx_sorted[i] = i;
    /* Simple selection sort for top local_k */
    for (int kk = 0; kk < local_k && kk < n_theta; kk++) {
        int min_idx = kk;
        for (int j = kk + 1; j < n_theta; j++) {
            if (vals[idx_sorted[j]] < vals[idx_sorted[min_idx]]) min_idx = j;
        }
        int tmp = idx_sorted[kk]; idx_sorted[kk] = idx_sorted[min_idx]; idx_sorted[min_idx] = tmp;
    }

    double phi = (sqrt(5.0) - 1.0) / 2.0;

    for (int ki = 0; ki < local_k && ki < n_theta; ki++) {
        int idx = idx_sorted[ki];
        double a = thetas[idx > 0 ? idx - 1 : 0];
        double b = thetas[idx < n_theta - 1 ? idx + 1 : n_theta - 1];
        if (b - a <= 1e-15) continue;

        double c = b - (b - a) * phi;
        double d = a + (b - a) * phi;
        Point rc[5], rd[5];
        rotate_points(points, rc, 5, c);
        rotate_points(points, rd, 5, d);
        double rc_min = INFINITY, rc_max = -INFINITY, rd_min = INFINITY, rd_max = -INFINITY;
        for (int j = 0; j < 5; j++) {
            if (rc[j].y < rc_min) rc_min = rc[j].y;
            if (rc[j].x > rc_max) rc_max = rc[j].x;
            if (rd[j].y < rd_min) rd_min = rd[j].y;
            if (rd[j].x > rd_max) rd_max = rd[j].x;
        }
        double fc = clearance_for_theta(rc, rc_min, rc_max, scale, eps_y, 5);
        double fd = clearance_for_theta(rd, rd_min, rd_max, scale, eps_y, 5);

        for (int iter = 0; iter < local_iters; iter++) {
            if (fc < fd) {
                b = d; d = c; fd = fc;
                c = b - (b - a) * phi;
                rotate_points(points, rc, 5, c);
                rc_min = INFINITY; rc_max = -INFINITY;
                for (int j = 0; j < 5; j++) { if (rc[j].y < rc_min) rc_min = rc[j].y; if (rc[j].x > rc_max) rc_max = rc[j].x; }
                fc = clearance_for_theta(rc, rc_min, rc_max, scale, eps_y, 5);
            } else {
                a = c; c = d; fc = fd;
                d = a + (b - a) * phi;
                rotate_points(points, rd, 5, d);
                rd_min = INFINITY; rd_max = -INFINITY;
                for (int j = 0; j < 5; j++) { if (rd[j].y < rd_min) rd_min = rd[j].y; if (rd[j].x > rd_max) rd_max = rd[j].x; }
                fd = clearance_for_theta(rd, rd_min, rd_max, scale, eps_y, 5);
            }
        }
        if (fc < best) best = fc;
        if (fd < best) best = fd;
    }

    free(vals);
    free(idx_sorted);
    return best;
}

static double max_scale(const Point *points, int n_theta, int bisection_iters, double eps_y) {
    double *thetas = malloc(n_theta * sizeof(double));
    RotData *precomp = malloc(n_theta * sizeof(RotData));
    for (int i = 0; i < n_theta; i++) {
        thetas[i] = (M_PI / 2.0) * i / (n_theta - 1);
        rotate_points(points, precomp[i].pts, 5, thetas[i]);
        precomp[i].miny = INFINITY; precomp[i].maxx = -INFINITY;
        for (int j = 0; j < 5; j++) {
            if (precomp[i].pts[j].y < precomp[i].miny) precomp[i].miny = precomp[i].pts[j].y;
            if (precomp[i].pts[j].x > precomp[i].maxx) precomp[i].maxx = precomp[i].pts[j].x;
        }
    }

    double lo = 0.0, hi = 2.0;
    /* Expand hi until infeasible */
    while (1) {
        double mc = min_clearance(points, precomp, thetas, n_theta, hi, eps_y, 3, 22);
        if (mc < -1e-13) break;
        hi *= 1.3;
        if (hi > 50.0) break;
    }

    for (int iter = 0; iter < bisection_iters; iter++) {
        double mid = (lo + hi) * 0.5;
        double mc = min_clearance(points, precomp, thetas, n_theta, mid, eps_y, 3, 22);
        if (mc >= -1e-13) lo = mid;
        else hi = mid;
    }

    free(thetas);
    free(precomp);
    return lo;
}

typedef struct { const char *mode; int n_theta; int bisection_iters; double eps_y; } Mode;

static double objective(double r, const Mode *m) {
    Point pts[5];
    if (!build_unit_pentagon(r, pts)) return -1.0;
    double s = max_scale(pts, m->n_theta, m->bisection_iters, m->eps_y);
    return base_area(r) * s * s;
}

static double golden_max(double (*f)(double, void*), void *ctx, double a, double b, int iters) {
    double phi = (sqrt(5.0) - 1.0) / 2.0;
    double c = b - (b - a) * phi;
    double d = a + (b - a) * phi;
    double fc = f(c, ctx);
    double fd = f(d, ctx);
    for (int i = 0; i < iters; i++) {
        if (fc > fd) {
            b = d; d = c; fd = fc;
            c = b - (b - a) * phi;
            fc = f(c, ctx);
        } else {
            a = c; c = d; fc = fd;
            d = a + (b - a) * phi;
            fd = f(d, ctx);
        }
    }
    return fc > fd ? c : d;
}

/* Cache for f_mid */
typedef struct {
    double *keys;
    double *vals;
    int count;
    int cap;
} Cache;

static Cache *cache_new(int cap) {
    Cache *c = malloc(sizeof(Cache));
    c->keys = malloc(cap * sizeof(double));
    c->vals = malloc(cap * sizeof(double));
    c->count = 0;
    c->cap = cap;
    return c;
}

static double cache_get(Cache *c, double key, double (*compute)(double, void*), void *ctx) {
    for (int i = 0; i < c->count; i++) {
        if (fabs(c->keys[i] - key) < 1e-15) return c->vals[i];
    }
    double val = compute(key, ctx);
    if (c->count < c->cap) {
        c->keys[c->count] = key;
        c->vals[c->count] = val;
        c->count++;
    }
    return val;
}

static Mode mid_mode = {"mid", 1400, 55, 1e-14};

static double f_mid_compute(double rr, void *ctx) {
    return objective(rr, &mid_mode);
}

static double f_mid_cached(double rr, void *ctx) {
    Cache *c = (Cache *)ctx;
    return cache_get(c, rr, f_mid_compute, NULL);
}

double p667_native(void) {
    /* 1) Coarse scan */
    double rmin = 0.75, rmax = 1.05;
    int steps = 240;
    Mode coarse_mode = {"coarse", 450, 35, 1e-12};
    double best_r = 0, best_val = -1.0;
    for (int i = 0; i <= steps; i++) {
        double r = rmin + (rmax - rmin) * i / steps;
        double val = objective(r, &coarse_mode);
        if (val > best_val) { best_val = val; best_r = r; }
    }

    /* 2) Golden refine with mid precision + cache */
    Cache *cache = cache_new(10000);
    double a = best_r - 0.02, b = best_r + 0.02;
    double r1 = golden_max(f_mid_cached, cache, a, b, 26);

    /* 3) Narrow refine */
    double a2 = r1 - 0.002, b2 = r1 + 0.002;
    double r2 = golden_max(f_mid_cached, cache, a2, b2, 35);

    /* 4) Final high precision */
    Mode fine_mode = {"fine", 8000, 75, 1e-15};
    double final_area = objective(r2, &fine_mode);

    free(cache->keys);
    free(cache->vals);
    free(cache);

    return final_area;
}
