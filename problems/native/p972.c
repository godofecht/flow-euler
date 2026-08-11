/* Project Euler 972 — Rational triangles on the unit sphere.
 *
 * Port of the reference Python solver.  The unit sphere x^2 + y^2 + z^2 = 1
 * is scaled by L = lcm(1..n) so that every rational point with denominator
 * <= n lands on an integer lattice.  Points are lifted to the cone
 * z = x^2 + y^2 + L^2 and triples are counted by collinearity of the
 * projection through the origin.
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

typedef int64_t i64;

static i64 gcd64(i64 a, i64 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { i64 t = a % b; a = b; b = t; }
    return a;
}

static i64 lcm64(i64 a, i64 b) {
    return a / gcd64(a, b) * b;
}

static int cmp_i64(const void *p, const void *q) {
    i64 a = *(const i64 *)p, b = *(const i64 *)q;
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}

typedef struct { i64 a, b; } Key;

static int cmp_key(const void *p, const void *q) {
    const Key *x = p, *y = q;
    if (x->a < y->a) return -1;
    if (x->a > y->a) return 1;
    if (x->b < y->b) return -1;
    if (x->b > y->b) return 1;
    return 0;
}

typedef struct { i64 x, y, z; } Point;

static i64      g_scale;
static i64     *g_coords;
static int      g_ncoords;
static Point   *g_points;
static int      g_npoints;

static void coordinate_values(int n) {
    g_scale = 1;
    for (int d = 1; d <= n; d++) g_scale = lcm64(g_scale, d);

    int cap = 256;
    g_coords = malloc((size_t)cap * sizeof(i64));
    g_ncoords = 0;
    g_coords[g_ncoords++] = 0;

    for (int d = 2; d <= n; d++) {
        i64 step = g_scale / d;
        for (int a = 1; a < d; a++) {
            if (gcd64((i64)a, (i64)d) == 1) {
                i64 v = (i64)a * step;
                if (g_ncoords + 2 > cap) {
                    cap *= 2;
                    g_coords = realloc(g_coords, (size_t)cap * sizeof(i64));
                }
                g_coords[g_ncoords++] = v;
                g_coords[g_ncoords++] = -v;
            }
        }
    }

    qsort(g_coords, (size_t)g_ncoords, sizeof(i64), cmp_i64);

    int w = 0;
    for (int i = 0; i < g_ncoords; i++) {
        if (i == 0 || g_coords[i] != g_coords[i - 1]) g_coords[w++] = g_coords[i];
    }
    g_ncoords = w;
}

static void build_lifted_points(int n) {
    coordinate_values(n);
    i64 scale2 = g_scale * g_scale;

    i64 *sq = malloc((size_t)g_ncoords * sizeof(i64));
    for (int i = 0; i < g_ncoords; i++) sq[i] = g_coords[i] * g_coords[i];

    int cap = 1 << 16;
    g_points = malloc((size_t)cap * sizeof(Point));
    g_npoints = 0;

    for (int xi = 0; xi < g_ncoords; xi++) {
        i64 x  = g_coords[xi];
        i64 x2 = sq[xi];
        for (int yi = 0; yi < g_ncoords; yi++) {
            i64 z0 = x2 + sq[yi];
            if (z0 < scale2) {
                if (g_npoints >= cap) {
                    cap *= 2;
                    g_points = realloc(g_points, (size_t)cap * sizeof(Point));
                }
                g_points[g_npoints].x = x;
                g_points[g_npoints].y = g_coords[yi];
                g_points[g_npoints].z = z0 + scale2;
                g_npoints++;
            }
        }
    }

    free(sq);
}

static i64 count_unordered_triples(void) {
    i64 total = 0;
    int m = g_npoints;

    Key *keys = malloc((size_t)m * sizeof(Key));

    for (int i = 0; i < m - 2; i++) {
        i64 px = g_points[i].x;
        i64 py = g_points[i].y;
        i64 pz = g_points[i].z;
        int nk = 0;

        if (px != 0) {
            for (int j = i + 1; j < m; j++) {
                i64 qx = g_points[j].x;
                i64 qy = g_points[j].y;
                i64 qz = g_points[j].z;
                i64 a = py * qx - px * qy;
                i64 b = pz * qx - px * qz;
                i64 g = gcd64(a, b);
                if (g != 0) { a /= g; b /= g; }
                if (a < 0 || (a == 0 && b < 0)) { a = -a; b = -b; }
                keys[nk].a = a;
                keys[nk].b = b;
                nk++;
            }
        } else {
            for (int j = i + 1; j < m; j++) {
                i64 qx = g_points[j].x;
                i64 qy = g_points[j].y;
                i64 qz = g_points[j].z;
                i64 a = qx;
                i64 b = pz * qy - py * qz;
                i64 g = gcd64(a, b);
                if (g != 0) { a /= g; b /= g; }
                if (a < 0 || (a == 0 && b < 0)) { a = -a; b = -b; }
                keys[nk].a = a;
                keys[nk].b = b;
                nk++;
            }
        }

        qsort(keys, (size_t)nk, sizeof(Key), cmp_key);

        i64 run_len = 1;
        Key prev = keys[0];
        for (int k = 1; k < nk; k++) {
            if (keys[k].a == prev.a && keys[k].b == prev.b) {
                run_len++;
            } else {
                if (run_len >= 2) total += run_len * (run_len - 1) / 2;
                prev = keys[k];
                run_len = 1;
            }
        }
        if (run_len >= 2) total += run_len * (run_len - 1) / 2;
    }

    free(keys);
    return total;
}

static i64 compute_T(int n) {
    build_lifted_points(n);
    i64 c = count_unordered_triples();
    i64 result = 6 * c;
    free(g_coords);
    free(g_points);
    return result;
}

long long p972_native(void) {
    return (long long)compute_T(12);
}
