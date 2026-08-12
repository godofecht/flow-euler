/* Project Euler 966: Triangle Circle Intersection
 *
 * sum I(a,b,c) for 1<=a<=b<=c<a+b, a+b+c<=200, where I(a,b,c) is the
 * largest intersection area between the triangle and a circle of equal
 * area, with the circle free to translate.
 *
 * Answer: 29337152.09
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef long long i64;

static double PI;
static const double EPS = 1e-12;

/* ---------- segment / circle intersection ---------- */

static int segment_circle_ts(double ax, double ay, double bx, double by,
                              double r2, double *t1_out, double *t2_out) {
    double dx = bx - ax, dy = by - ay;
    double qa = dx * dx + dy * dy;
    if (qa < 1e-18) { *t1_out = *t2_out = 0; return 0; }

    double qb = 2.0 * (ax * dx + ay * dy);
    double qc = ax * ax + ay * ay - r2;
    double disc = qb * qb - 4.0 * qa * qc;
    if (disc < -1e-12) { *t1_out = *t2_out = 0; return 0; }
    if (disc < 0) disc = 0;

    double sdisc = sqrt(disc);
    double inv2a = 0.5 / qa;
    double t1 = (-qb - sdisc) * inv2a;
    double t2 = (-qb + sdisc) * inv2a;

    double ts[2];
    int nts = 0;
    if (-EPS <= t1 && t1 <= 1.0 + EPS) {
        ts[nts] = (t1 < 0) ? 0 : (t1 > 1 ? 1 : t1);
        nts++;
    }
    if (-EPS <= t2 && t2 <= 1.0 + EPS) {
        double t2c = (t2 < 0) ? 0 : (t2 > 1 ? 1 : t2);
        if (nts == 0 || fabs(t2c - ts[0]) > 1e-11)
            ts[nts++] = t2c;
    }

    if (nts == 0) { *t1_out = *t2_out = 0; return 0; }
    if (nts == 1) { *t1_out = ts[0]; *t2_out = 0; return 1; }
    if (ts[0] > ts[1]) { double tmp = ts[0]; ts[0] = ts[1]; ts[1] = tmp; }
    *t1_out = ts[0]; *t2_out = ts[1];
    return 2;
}

static double tri_or_sector(double px, double py, double qx, double qy, double r2) {
    double cross = px * qy - py * qx;
    double p2 = px * px + py * py;
    double q2 = qx * qx + qy * qy;
    if (p2 <= r2 + 1e-12 && q2 <= r2 + 1e-12)
        return 0.5 * cross;
    double dot = px * qx + py * qy;
    double ang = atan2(cross, dot);
    return 0.5 * r2 * ang;
}

static void project_to_circle(double x, double y, double r,
                               double *ox, double *oy) {
    double d = hypot(x, y);
    if (d < 1e-18) { *ox = x; *oy = y; return; }
    double s = r / d;
    *ox = x * s; *oy = y * s;
}

static double edge_contrib(double ax, double ay, double bx, double by, double r2) {
    double t1, t2;
    int k = segment_circle_ts(ax, ay, bx, by, r2, &t1, &t2);
    if (k == 0)
        return tri_or_sector(ax, ay, bx, by, r2);

    double dx = bx - ax, dy = by - ay;
    double r = sqrt(r2);

    if (k == 1) {
        double ix = ax + dx * t1, iy = ay + dy * t1;
        project_to_circle(ix, iy, r, &ix, &iy);
        return tri_or_sector(ax, ay, ix, iy, r2)
             + tri_or_sector(ix, iy, bx, by, r2);
    }

    double i1x = ax + dx * t1, i1y = ay + dy * t1;
    double i2x = ax + dx * t2, i2y = ay + dy * t2;
    project_to_circle(i1x, i1y, r, &i1x, &i1y);
    project_to_circle(i2x, i2y, r, &i2x, &i2y);
    return tri_or_sector(ax, ay, i1x, i1y, r2)
         + tri_or_sector(i1x, i1y, i2x, i2y, r2)
         + tri_or_sector(i2x, i2y, bx, by, r2);
}

static double tri_circle_area(double ax, double ay, double bx, double by,
                               double cx, double cy,
                               double ox, double oy, double r) {
    double r2 = r * r;
    double a1x = ax - ox, a1y = ay - oy;
    double b1x = bx - ox, b1y = by - oy;
    double c1x = cx - ox, c1y = cy - oy;
    double total = 0;
    total += edge_contrib(a1x, a1y, b1x, b1y, r2);
    total += edge_contrib(b1x, b1y, c1x, c1y, r2);
    total += edge_contrib(c1x, c1y, a1x, a1y, r2);
    return fabs(total);
}

/* ---------- triangle construction ---------- */

static void tri_coords(int a, int b, int c,
                        double *ax, double *ay, double *bx, double *by,
                        double *cx, double *cy) {
    *ax = 0; *ay = 0;
    *bx = (double)c; *by = 0;
    double x = (double)(b * b + c * c - a * a) / (2.0 * c);
    double y2 = (double)(b * b) - x * x;
    if (y2 < 0 && y2 > -1e-12) y2 = 0;
    *cy = (y2 > 0) ? sqrt(y2) : 0;
    *cx = x;
}

static double tri_area(double ax, double ay, double bx, double by,
                        double cx, double cy) {
    return fabs(0.5 * ((bx - ax) * (cy - ay) - (by - ay) * (cx - ax)));
}

static int circumcenter(double ax, double ay, double bx, double by,
                         double cx, double cy, double *ux, double *uy) {
    double d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    if (fabs(d) < 1e-15) return 0;
    double a2 = ax * ax + ay * ay;
    double b2 = bx * bx + by * by;
    double c2 = cx * cx + cy * cy;
    *ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d;
    *uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d;
    return 1;
}

/* ---------- maximization (pattern search) ---------- */

static double maximize_intersection(int a, int b, int c) {
    double ax, ay, bx, by, cx, cy;
    tri_coords(a, b, c, &ax, &ay, &bx, &by, &cx, &cy);
    double area = tri_area(ax, ay, bx, by, cx, cy);
    if (area <= 0) return 0;

    double r = sqrt(area / PI);

    double xmin = fmin(ax, fmin(bx, cx)) - r;
    double xmax = fmax(ax, fmax(bx, cx)) + r;
    double ymin = fmin(ay, fmin(by, cy)) - r;
    double ymax = fmax(ay, fmax(by, cy)) + r;

    double cenx = (ax + bx + cx) / 3.0;
    double ceny = (ay + by + cy) / 3.0;

    int per = a + b + c;
    double incx = (a * ax + b * bx + c * cx) / (double)per;
    double incy = (a * ay + b * by + c * cy) / (double)per;

    double ccx, ccy;
    int has_cc = circumcenter(ax, ay, bx, by, cx, cy, &ccx, &ccy);

    double starts[16][2];
    int ns = 0;
    starts[ns][0] = cenx;          starts[ns][1] = ceny;          ns++;
    starts[ns][0] = incx;          starts[ns][1] = incy;          ns++;
    starts[ns][0] = ax;            starts[ns][1] = ay;            ns++;
    starts[ns][0] = bx;            starts[ns][1] = by;            ns++;
    starts[ns][0] = cx;            starts[ns][1] = cy;            ns++;
    starts[ns][0] = (ax + bx)*0.5; starts[ns][1] = (ay + by)*0.5; ns++;
    starts[ns][0] = (bx + cx)*0.5; starts[ns][1] = (by + cy)*0.5; ns++;
    starts[ns][0] = (cx + ax)*0.5; starts[ns][1] = (cy + ay)*0.5; ns++;
    if (has_cc) { starts[ns][0] = ccx; starts[ns][1] = ccy; ns++; }

    /* clamp helper */
    double bestx = cenx, besty = ceny;
    if (bestx < xmin) bestx = xmin; else if (bestx > xmax) bestx = xmax;
    if (besty < ymin) besty = ymin; else if (besty > ymax) besty = ymax;

    double bestv = tri_circle_area(ax, ay, bx, by, cx, cy, bestx, besty, r);
    if (bestv < 0) bestv = 0;
    if (bestv > area) bestv = area;

    for (int s = 0; s < ns; s++) {
        double sx = starts[s][0], sy = starts[s][1];
        if (sx < xmin) sx = xmin; else if (sx > xmax) sx = xmax;
        if (sy < ymin) sy = ymin; else if (sy > ymax) sy = ymax;
        double v = tri_circle_area(ax, ay, bx, by, cx, cy, sx, sy, r);
        if (v < 0) v = 0;
        if (v > area) v = area;
        if (v > bestv) { bestv = v; bestx = sx; besty = sy; }
    }

    double span = fmax(xmax - xmin, ymax - ymin);
    double step = fmax(span, r);

    static const double dirs[8][2] = {
        {1,0},{-1,0},{0,1},{0,-1},
        {1,1},{1,-1},{-1,1},{-1,-1}
    };

    double tol_step = fmax(1e-7, 1e-7 * fmax(1.0, r));

    while (step > tol_step) {
        int moved_any = 0;
        while (1) {
            int improved = 0;
            double bx0 = bestx, by0 = besty, bv0 = bestv;
            for (int d = 0; d < 8; d++) {
                double px = bestx + dirs[d][0] * step;
                double py = besty + dirs[d][1] * step;
                if (px < xmin) px = xmin; else if (px > xmax) px = xmax;
                if (py < ymin) py = ymin; else if (py > ymax) py = ymax;
                double v = tri_circle_area(ax, ay, bx, by, cx, cy, px, py, r);
                if (v < 0) v = 0;
                if (v > area) v = area;
                if (v > bv0 + 1e-13) {
                    bv0 = v; bx0 = px; by0 = py;
                    improved = 1;
                }
            }
            if (improved) {
                bestv = bv0; bestx = bx0; besty = by0;
                moved_any = 1;
            } else {
                break;
            }
        }
        if (!moved_any) step *= 0.5;
    }

    return bestv;
}

/* ---------- shape coefficients ---------- */

typedef struct { int a, b, c; i64 g2; } TriTuple;

static int cmp_trituple(const void *p, const void *q) {
    const TriTuple *x = p, *y = q;
    if (x->a != y->a) return x->a - y->a;
    if (x->b != y->b) return x->b - y->b;
    return x->c - y->c;
}

static int gcd_int(int a, int b) {
    while (b) { int t = a % b; a = b; b = t; }
    return a;
}

double p966_native(void) {
    PI = acos(-1.0);
    int limit = 200;

    /* Collect (primitive_shape, g^2) tuples for every triangle. */
    int cap = 300000;
    TriTuple *tuples = malloc((size_t)cap * sizeof(TriTuple));
    int nt = 0;

    for (int a = 1; a <= limit; a++) {
        for (int b = a; b <= limit; b++) {
            int maxc = (a + b - 1 < limit - a - b) ? a + b - 1 : limit - a - b;
            if (maxc < b) continue;
            for (int c = b; c <= maxc; c++) {
                int g = gcd_int(gcd_int(a, b), c);
                tuples[nt].a = a / g;
                tuples[nt].b = b / g;
                tuples[nt].c = c / g;
                tuples[nt].g2 = (i64)g * g;
                nt++;
            }
        }
    }

    qsort(tuples, (size_t)nt, sizeof(TriTuple), cmp_trituple);

    /* Merge duplicates and accumulate weights with Kahan summation. */
    double total = 0.0, corr = 0.0;
    int i = 0;
    while (i < nt) {
        int a = tuples[i].a, b = tuples[i].b, c = tuples[i].c;
        i64 weight = 0;
        while (i < nt && tuples[i].a == a &&
               tuples[i].b == b && tuples[i].c == c) {
            weight += tuples[i].g2;
            i++;
        }
        double Ival = maximize_intersection(a, b, c);
        double term = Ival * (double)weight;
        double y = term - corr;
        double t = total + y;
        corr = (t - total) - y;
        total = t;
    }

    free(tuples);
    return total;
}
