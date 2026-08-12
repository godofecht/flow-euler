// Project Euler 842: Irregular Star Polygons
// Port of Python reference solver to C.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MOD 1000000007LL

static long long modpow(long long b, long long e, long long m) {
    long long r = 1; b %= m;
    while (e > 0) { if (e&1) r = r*b%m; b = b*b%m; e >>= 1; }
    return r;
}
static long long modinv(long long a, long long m) { return modpow(a, m-2, m); }

static long long fact_mod[61], inv_fact_mod[61];
static void init_factorials() {
    fact_mod[0] = 1;
    for (int i = 1; i <= 60; i++) fact_mod[i] = fact_mod[i-1] * i % MOD;
    for (int i = 0; i <= 60; i++) inv_fact_mod[i] = modinv(fact_mod[i], MOD);
}
static long long comb_mod(int n, int k) {
    if (k < 0 || k > n) return 0;
    return fact_mod[n] * inv_fact_mod[k] % MOD * inv_fact_mod[n-k] % MOD;
}

static long long cycles_at_least_two_edges(int n, int m) {
    if (m < 2) return 0;
    long long inv2 = modinv(2, MOD);
    long long total = fact_mod[n-1] * inv2 % MOD;
    long long c0 = 0;
    for (int k = 0; k <= m; k++) {
        long long ak;
        if (k == 0) ak = total;
        else ak = modpow(2, k-1, MOD) * fact_mod[n-k-1] % MOD;
        long long term = comb_mod(m, k) * ak % MOD;
        if (k % 2 == 0) c0 = (c0 + term) % MOD;
        else c0 = (c0 - term + MOD) % MOD;
    }
    long long c1 = 0;
    for (int k = 1; k <= m; k++) {
        long long ak = modpow(2, k-1, MOD) * fact_mod[n-k-1] % MOD;
        long long term = (long long)k * comb_mod(m, k) % MOD * ak % MOD;
        if ((k-1) % 2 == 0) c1 = (c1 + term) % MOD;
        else c1 = (c1 - term + MOD) % MOD;
    }
    return (total - c0 - c1 % MOD + 2*MOD) % MOD;
}

static long long isqrt_ll(long long n) {
    if (n < 0) return -1; if (n == 0) return 0;
    long long x = (long long)sqrt((double)n);
    while (x*x > n) x--; while ((x+1)*(x+1) <= n) x++; return x;
}

static int inv_triangular(int q) {
    long long disc = 1 + 8LL*q; long long r = isqrt_ll(disc);
    if (r*r != disc) return -1; if ((1+r)%2 != 0) return -1; return (int)((1+r)/2);
}

// ---- Fast path: quantized hash ----
#define SCALE 100000000000LL
#define MAX_POINTS 600000
#define HASH_SIZE (1 << 20)
#define HASH_MASK (HASH_SIZE - 1)

typedef struct { long long qx, qy; } PointKey;
typedef struct { PointKey key; int count; } PointEntry;

static PointEntry points[MAX_POINTS];
static int num_points;
static int hash_buckets[HASH_SIZE];
static int hash_next[MAX_POINTS];

static unsigned int hash_key(long long qx, long long qy) {
    unsigned long long h = (unsigned long long)qx * 1000000007ULL + (unsigned long long)qy;
    h ^= h >> 32;
    return (unsigned int)(h & HASH_MASK);
}
static int find_point(long long qx, long long qy) {
    unsigned int h = hash_key(qx, qy);
    int idx = hash_buckets[h];
    while (idx != -1) {
        if (points[idx].key.qx == qx && points[idx].key.qy == qy) return idx;
        idx = hash_next[idx];
    }
    return -1;
}
static void add_point_fast(long long qx, long long qy) {
    int idx = find_point(qx, qy);
    if (idx >= 0) { points[idx].count++; return; }
    idx = num_points++;
    points[idx].key.qx = qx; points[idx].key.qy = qy; points[idx].count = 1;
    unsigned int h = hash_key(qx, qy);
    hash_next[idx] = hash_buckets[h]; hash_buckets[h] = idx;
}

// ---- Slow fallback: spatial clustering ----
#define CELL 1e-6
#define EPS2 1e-18
#define MAX_REPS 600000

static double reps_x[MAX_REPS], reps_y[MAX_REPS];
static int reps_cnt[MAX_REPS];
static int num_reps;
static int grid_buckets[HASH_SIZE];
static int grid_next[MAX_REPS];

static unsigned int grid_key(double x, double y) {
    int ix = (int)floor(x / CELL);
    int iy = (int)floor(y / CELL);
    unsigned long long h = (unsigned long long)ix * 1000000007ULL + (unsigned long long)iy;
    h ^= h >> 32;
    return (unsigned int)(h & HASH_MASK);
}

static void add_point_slow(double x, double y) {
    int ix = (int)floor(x / CELL);
    int iy = (int)floor(y / CELL);
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            int cx = ix + dx, cy = iy + dy;
            unsigned long long h2 = (unsigned long long)cx * 1000000007ULL + (unsigned long long)cy;
            h2 ^= h2 >> 32;
            unsigned int h = (unsigned int)(h2 & HASH_MASK);
            int idx = grid_buckets[h];
            while (idx != -1) {
                double ddx = x - reps_x[idx];
                double ddy = y - reps_y[idx];
                if (ddx*ddx + ddy*ddy <= EPS2) {
                    reps_cnt[idx]++;
                    return;
                }
                idx = grid_next[idx];
            }
        }
    }
    int idx = num_reps++;
    reps_x[idx] = x; reps_y[idx] = y; reps_cnt[idx] = 1;
    unsigned int h = grid_key(x, y);
    grid_next[idx] = grid_buckets[h]; grid_buckets[h] = idx;
}

// Compute intersection points for given n
// Returns 1 if fast path succeeded (all counts triangular), 0 if fallback needed
static int compute_points_fast(int n, double *xs, double *ys) {
    num_points = 0;
    for (int i = 0; i < HASH_SIZE; i++) hash_buckets[i] = -1;

    for (int a = 0; a < n-3; a++) {
        double x1=xs[a], y1=ys[a];
        for (int b = a+1; b < n-2; b++) {
            double x3=xs[b], y3=ys[b];
            double abx=x3-x1, aby=y3-y1;
            for (int c = b+1; c < n-1; c++) {
                double x2=xs[c], y2=ys[c];
                double dx12=x2-x1, dy12=y2-y1;
                for (int d = c+1; d < n; d++) {
                    double x4=xs[d], y4=ys[d];
                    double dx34=x4-x3, dy34=y4-y3;
                    double denom=dx12*dy34-dy12*dx34;
                    if (fabs(denom)<1e-18) continue;
                    double t=(abx*dy34-aby*dx34)/denom;
                    double px=x1+t*dx12, py=y1+t*dy12;
                    long long qx=(long long)llround(px*SCALE), qy=(long long)llround(py*SCALE);
                    add_point_fast(qx, qy);
                }
            }
        }
    }
    // Check all counts are triangular
    for (int i = 0; i < num_points; i++) {
        if (inv_triangular(points[i].count) < 0) return 0;
    }
    return 1;
}

static void compute_points_slow(int n, double *xs, double *ys) {
    num_reps = 0;
    for (int i = 0; i < HASH_SIZE; i++) grid_buckets[i] = -1;

    for (int a = 0; a < n-3; a++) {
        double x1=xs[a], y1=ys[a];
        for (int b = a+1; b < n-2; b++) {
            double x3=xs[b], y3=ys[b];
            double abx=x3-x1, aby=y3-y1;
            for (int c = b+1; c < n-1; c++) {
                double x2=xs[c], y2=ys[c];
                double dx12=x2-x1, dy12=y2-y1;
                for (int d = c+1; d < n; d++) {
                    double x4=xs[d], y4=ys[d];
                    double dx34=x4-x3, dy34=y4-y3;
                    double denom=dx12*dy34-dy12*dx34;
                    if (fabs(denom)<1e-18) continue;
                    double t=(abx*dy34-aby*dx34)/denom;
                    double px=x1+t*dx12, py=y1+t*dy12;
                    add_point_slow(px, py);
                }
            }
        }
    }
}

static long long T_mod(int n) {
    if (n < 4) return 0;

    double tau = 2.0 * M_PI;
    double xs[60], ys[60];
    for (int k = 0; k < n; k++) { xs[k] = cos(tau*k/n); ys[k] = sin(tau*k/n); }

    int dist_m[100], dist_cnt[100], dist_size = 0;
    memset(dist_m, 0, sizeof(dist_m)); memset(dist_cnt, 0, sizeof(dist_cnt));

    if (compute_points_fast(n, xs, ys)) {
        for (int i = 0; i < num_points; i++) {
            int m = inv_triangular(points[i].count);
            int found = -1;
            for (int j = 0; j < dist_size; j++) if (dist_m[j] == m) { found = j; break; }
            if (found < 0) { dist_m[dist_size] = m; dist_cnt[dist_size] = 1; dist_size++; }
            else dist_cnt[found]++;
        }
    } else {
        compute_points_slow(n, xs, ys);
        for (int i = 0; i < num_reps; i++) {
            int m = inv_triangular(reps_cnt[i]);
            if (m < 0) continue;
            int found = -1;
            for (int j = 0; j < dist_size; j++) if (dist_m[j] == m) { found = j; break; }
            if (found < 0) { dist_m[dist_size] = m; dist_cnt[dist_size] = 1; dist_size++; }
            else dist_cnt[found]++;
        }
    }

    long long total = 0;
    for (int i = 0; i < dist_size; i++) {
        long long g = cycles_at_least_two_edges(n, dist_m[i]);
        total = (total + (long long)dist_cnt[i] * g) % MOD;
    }
    return total;
}

long long p842_native(void) {
    init_factorials();
    long long ans = 0;
    for (int n = 3; n <= 60; n++) {
        ans = (ans + T_mod(n)) % MOD;
    }
    return ans;
}
