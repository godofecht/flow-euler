#include <stdint.h>
#include <string.h>
#include <stdlib.h>

typedef long long i64;
typedef __int128 i128;

#define MOD 1000000007LL
#define M 8
#define W 64
#define SUM_OFF 2000
#define SUM_RNG 4001

static i64 powmod(i64 a, i64 b, i64 m) {
    i64 r = 1; a %= m; if (a < 0) a += m;
    while (b > 0) {
        if (b & 1) r = (i128)r * a % m;
        a = (i128)a * a % m;
        b >>= 1;
    }
    return r;
}

/* ---------- reduced_hook ---------- */

static void reduced_hook(i64 a, i64 b, i64 k, i64 *M_out, i64 *N_out) {
    i64 d = 0;
    for (i64 j = 0; j < k; j++) {
        i64 row_len = (k - j) * b;
        i64 start = j * a + 1;
        if (row_len < start) continue;
        i64 end = (j + 1) * a;
        i64 cand = (end < row_len) ? end : row_len;
        if (cand > d) d = cand;
    }
    i64 block_of_row_d = (d - 1) / a;
    i64 lambda_d = (k - block_of_row_d) * b;
    i64 Mv = lambda_d - d + 1;
    i64 need_blocks = (d + b - 1) / b;
    i64 last_block = k - need_blocks;
    i64 col_height = (last_block + 1) * a;
    i64 Nv = col_height - d + 1;
    *M_out = Mv;
    *N_out = Nv;
}

/* ---------- collect counts ---------- */

typedef struct { i64 t; i64 R; i64 count; } hot_type_t;

static hot_type_t hot_types[2000];
static int num_hot_types;

static i64 int_counts[2001];  /* int_counts[v + 1000] */
static i64 int_vals[2001];
static int num_int_vals;

static int cmp_hot(const void *a, const void *b) {
    const hot_type_t *ha = (const hot_type_t *)a;
    const hot_type_t *hb = (const hot_type_t *)b;
    if (ha->t != hb->t) return (ha->t < hb->t) ? 1 : -1;
    if (ha->R != hb->R) return (ha->R < hb->R) ? -1 : 1;
    return 0;
}

static void collect_counts(void) {
    memset(int_counts, 0, sizeof(int_counts));
    static hot_type_t temp_hot[50000];
    int temp_count = 0;

    for (i64 a = 1; a < W - 1; a++) {
        for (i64 b = 1; b < W - a; b++) {
            i64 max_k = W - a - b;
            if (max_k < 1) continue;
            for (i64 k = 1; k <= max_k; k++) {
                i64 Mv, Nv;
                reduced_hook(a, b, k, &Mv, &Nv);
                if (Nv == 1) {
                    i64 v = Mv - 1;
                    int_counts[v + 1000]++;
                } else if (Mv == 1) {
                    i64 v = -(Nv - 1);
                    int_counts[v + 1000]++;
                } else {
                    i64 t = Mv + Nv - 4;
                    i64 R = -(Nv - 2);
                    temp_hot[temp_count].t = t;
                    temp_hot[temp_count].R = R;
                    temp_hot[temp_count].count = 1;
                    temp_count++;
                }
            }
        }
    }

    qsort(temp_hot, temp_count, sizeof(hot_type_t), cmp_hot);

    num_hot_types = 0;
    for (int i = 0; i < temp_count; i++) {
        if (num_hot_types > 0 &&
            hot_types[num_hot_types - 1].t == temp_hot[i].t &&
            hot_types[num_hot_types - 1].R == temp_hot[i].R) {
            hot_types[num_hot_types - 1].count++;
        } else {
            hot_types[num_hot_types] = temp_hot[i];
            num_hot_types++;
        }
    }

    num_int_vals = 0;
    for (int v = -1000; v <= 1000; v++) {
        if (int_counts[v + 1000] > 0) {
            int_vals[num_int_vals] = v;
            num_int_vals++;
        }
    }
}

/* ---------- DP ---------- */

static i64 hot_a[9][2][SUM_RNG], hot_b[9][2][SUM_RNG];
static i64 int_a[9][SUM_RNG], int_b[9][SUM_RNG];

static i64 hot_min_a[9][2], hot_max_a[9][2];
static i64 hot_min_b[9][2], hot_max_b[9][2];
static i64 int_min_a[9], int_max_a[9];
static i64 int_min_b[9], int_max_b[9];

long long p923_native(void) {
    i64 fact[M + 1], invfact[M + 1];
    fact[0] = 1;
    for (int i = 1; i <= M; i++) fact[i] = fact[i - 1] * i % MOD;
    invfact[M] = powmod(fact[M], MOD - 2, MOD);
    for (int i = M; i >= 1; i--) invfact[i - 1] = invfact[i] * i % MOD;

    collect_counts();

    /* Hot DP */
    memset(hot_a, 0, sizeof(hot_a));
    for (int u = 0; u <= M; u++) for (int p = 0; p < 2; p++) {
        hot_min_a[u][p] = SUM_RNG; hot_max_a[u][p] = -1;
    }
    hot_a[0][0][SUM_OFF] = 1;
    hot_min_a[0][0] = SUM_OFF; hot_max_a[0][0] = SUM_OFF;

    for (int ht = 0; ht < num_hot_types; ht++) {
        i64 t = hot_types[ht].t;
        i64 R = hot_types[ht].R;
        i64 c = hot_types[ht].count;

        i64 poly[M + 1];
        poly[0] = 1;
        { i64 p = 1;
          for (int k = 1; k <= M; k++) { p = p * c % MOD; poly[k] = p * invfact[k] % MOD; }
        }

        memset(hot_b, 0, sizeof(hot_b));
        for (int u = 0; u <= M; u++) for (int p = 0; p < 2; p++) {
            hot_min_b[u][p] = SUM_RNG; hot_max_b[u][p] = -1;
        }

        for (int used = 0; used <= M; used++) {
            for (int parity = 0; parity < 2; parity++) {
                if (hot_max_a[used][parity] < 0) continue;
                for (int idx = hot_min_a[used][parity]; idx <= hot_max_a[used][parity]; idx++) {
                    i64 coeff = hot_a[used][parity][idx];
                    if (coeff == 0) continue;
                    for (int k = 0; k <= M - used; k++) {
                        i64 mult = poly[k];
                        if (mult == 0) continue;
                        int right_turns = (k + 1 - parity) / 2;
                        i64 delta = (i64)k * R + (i64)right_turns * t;
                        int nu = used + k;
                        int np = parity ^ (k & 1);
                        int nidx = idx + (int)delta;
                        if (nidx < 0 || nidx >= SUM_RNG) continue;
                        hot_b[nu][np][nidx] = (hot_b[nu][np][nidx] + coeff * mult) % MOD;
                        if (nidx < hot_min_b[nu][np]) hot_min_b[nu][np] = nidx;
                        if (nidx > hot_max_b[nu][np]) hot_max_b[nu][np] = nidx;
                    }
                }
            }
        }

        memcpy(hot_a, hot_b, sizeof(hot_b));
        memcpy(hot_min_a, hot_min_b, sizeof(hot_min_b));
        memcpy(hot_max_a, hot_max_b, sizeof(hot_max_b));
    }

    /* Int DP */
    memset(int_a, 0, sizeof(int_a));
    for (int u = 0; u <= M; u++) { int_min_a[u] = SUM_RNG; int_max_a[u] = -1; }
    int_a[0][SUM_OFF] = 1;
    int_min_a[0] = SUM_OFF; int_max_a[0] = SUM_OFF;

    for (int iv = 0; iv < num_int_vals; iv++) {
        i64 v = int_vals[iv];
        i64 c = int_counts[v + 1000];

        i64 poly[M + 1];
        poly[0] = 1;
        { i64 p = 1;
          for (int k = 1; k <= M; k++) { p = p * c % MOD; poly[k] = p * invfact[k] % MOD; }
        }

        memset(int_b, 0, sizeof(int_b));
        for (int u = 0; u <= M; u++) { int_min_b[u] = SUM_RNG; int_max_b[u] = -1; }

        for (int used = 0; used <= M; used++) {
            if (int_max_a[used] < 0) continue;
            for (int idx = int_min_a[used]; idx <= int_max_a[used]; idx++) {
                i64 coeff = int_a[used][idx];
                if (coeff == 0) continue;
                for (int k = 0; k <= M - used; k++) {
                    i64 mult = poly[k];
                    if (mult == 0) continue;
                    int nu = used + k;
                    int nidx = idx + (int)(k * v);
                    if (nidx < 0 || nidx >= SUM_RNG) continue;
                    int_b[nu][nidx] = (int_b[nu][nidx] + coeff * mult) % MOD;
                    if (nidx < int_min_b[nu]) int_min_b[nu] = nidx;
                    if (nidx > int_max_b[nu]) int_max_b[nu] = nidx;
                }
            }
        }

        memcpy(int_a, int_b, sizeof(int_b));
        memcpy(int_min_a, int_min_b, sizeof(int_min_b));
        memcpy(int_max_a, int_max_b, sizeof(int_max_b));
    }

    /* Combine */
    i64 multiset_count = 0;
    for (int j = 0; j <= M; j++) {
        for (int parity = 0; parity < 2; parity++) {
            if (hot_max_a[j][parity] < 0) continue;
            if (int_max_a[M - j] < 0) continue;
            for (int hidx = hot_min_a[j][parity]; hidx <= hot_max_a[j][parity]; hidx++) {
                i64 ch = hot_a[j][parity][hidx];
                if (ch == 0) continue;
                for (int iidx = int_min_a[M - j]; iidx <= int_max_a[M - j]; iidx++) {
                    i64 ci = int_a[M - j][iidx];
                    if (ci == 0) continue;
                    i64 total = (i64)(hidx - SUM_OFF) + (i64)(iidx - SUM_OFF);
                    if (total > 0 || (total == 0 && parity == 1)) {
                        multiset_count = (multiset_count + ch * ci) % MOD;
                    }
                }
            }
        }
    }

    return multiset_count * fact[M] % MOD;
}
