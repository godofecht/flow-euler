/* Project Euler 818: SET
   F(12) = sum over 12-card subsets C of S(C)^4
   where S(C) is the number of SETs (affine lines) in C.
   Port of the Python reference solver. */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef __int128 i128;

/* 81-bit masks using __int128 */
typedef i128 mask_t;

static int popcount128(mask_t x) {
    return __builtin_popcountll((unsigned long long)x) +
           __builtin_popcountll((unsigned long long)(x >> 64));
}

/* Binomial coefficient C(n, k) */
static long long nCk(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    long long r = 1;
    for (int i = 0; i < k; i++) {
        r = r * (n - i) / (i + 1);
    }
    return r;
}

/* Base-3 digits for 81 points */
static int d0[81], d1[81], d2[81], d3[81];

static void init_digits(void) {
    for (int pid = 0; pid < 81; pid++) {
        int x = pid;
        d0[pid] = x % 3; x /= 3;
        d1[pid] = x % 3; x /= 3;
        d2[pid] = x % 3; x /= 3;
        d3[pid] = x % 3;
    }
}

static int third_point(int i, int j) {
    int a0 = (d0[i] + d0[j]) % 3;
    int a1 = (d1[i] + d1[j]) % 3;
    int a2 = (d2[i] + d2[j]) % 3;
    int a3 = (d3[i] + d3[j]) % 3;
    int c0 = (3 - a0) % 3;
    int c1 = (3 - a1) % 3;
    int c2 = (3 - a2) % 3;
    int c3 = (3 - a3) % 3;
    return c0 + 3 * c1 + 9 * c2 + 27 * c3;
}

static int vec_key(int v0, int v1, int v2, int v3) {
    return v0 + 3 * v1 + 9 * v2 + 27 * v3;
}

static int neg_key[81];

static void init_neg_keys(void) {
    for (int key = 0; key < 81; key++) {
        int v0 = (2 * (key % 3)) % 3;
        int v1 = (2 * ((key / 3) % 3)) % 3;
        int v2 = (2 * ((key / 9) % 3)) % 3;
        int v3 = (2 * ((key / 27) % 3)) % 3;
        neg_key[key] = vec_key(v0, v1, v2, v3);
    }
}

/* Line triples: sorted (a, b, c) with a < b < c */
typedef struct { int a, b, c; } Triple;

static int cmp_triple(const void *pa, const void *pb) {
    const Triple *a = (const Triple*)pa;
    const Triple *b = (const Triple*)pb;
    if (a->a != b->a) return a->a - b->a;
    if (a->b != b->b) return a->b - b->b;
    return a->c - b->c;
}

static mask_t line_masks[1080];
static int line_dirs[1080];

static void build_geometry(void) {
    init_digits();
    init_neg_keys();

    /* Direction index */
    int dir_index[81];
    memset(dir_index, -1, sizeof(dir_index));
    int dirs[40];
    int ndirs = 0;

    for (int key = 1; key < 81; key++) {
        int canon = (key < neg_key[key]) ? key : neg_key[key];
        if (dir_index[canon] == -1) {
            dir_index[canon] = ndirs;
            dirs[ndirs] = canon;
            ndirs++;
        }
    }

    /* Enumerate all lines from all unordered point pairs */
    Triple *triples = (Triple*)malloc(12000 * sizeof(Triple));
    int ntriples = 0;

    for (int i = 0; i < 81; i++) {
        for (int j = i + 1; j < 81; j++) {
            int k = third_point(i, j);
            if (k == i || k == j) continue;
            int a, b, c;
            /* sort (i, j, k) */
            a = i; b = j; c = k;
            if (a > b) { int t = a; a = b; b = t; }
            if (b > c) { int t = b; b = c; c = t; }
            if (a > b) { int t = a; a = b; b = t; }
            triples[ntriples].a = a;
            triples[ntriples].b = b;
            triples[ntriples].c = c;
            ntriples++;
        }
    }

    /* Sort and deduplicate */
    qsort(triples, ntriples, sizeof(Triple), cmp_triple);
    int unique = 0;
    for (int i = 0; i < ntriples; i++) {
        if (i == 0 || triples[i].a != triples[i-1].a || triples[i].b != triples[i-1].b || triples[i].c != triples[i-1].c) {
            triples[unique++] = triples[i];
        }
    }
    ntriples = unique; /* should be 1080 */

    for (int idx = 0; idx < ntriples; idx++) {
        int a = triples[idx].a, b = triples[idx].b, c = triples[idx].c;
        line_masks[idx] = ((mask_t)1 << a) | ((mask_t)1 << b) | ((mask_t)1 << c);

        /* Direction from a -> b */
        int v0 = (d0[b] - d0[a] + 3) % 3;
        int v1 = (d1[b] - d1[a] + 3) % 3;
        int v2 = (d2[b] - d2[a] + 3) % 3;
        int v3 = (d3[b] - d3[a] + 3) % 3;
        int key = vec_key(v0, v1, v2, v3);
        if (key == 0) {
            v0 = (d0[c] - d0[a] + 3) % 3;
            v1 = (d1[c] - d1[a] + 3) % 3;
            v2 = (d2[c] - d2[a] + 3) % 3;
            v3 = (d3[c] - d3[a] + 3) % 3;
            key = vec_key(v0, v1, v2, v3);
        }
        int canon = (key < neg_key[key]) ? key : neg_key[key];
        line_dirs[idx] = dir_index[canon];
    }

    free(triples);
}

static void compute_Ak(i128 *A_k) {
    int L = 1080;

    /* Classes of ordered line pairs (i,j):
       0: i==j (same line), union size 3
       1: intersecting distinct lines, union size 5
       2: disjoint parallel distinct lines, union size 6
       3: disjoint non-parallel (skew), union size 6 */

    /* We need to store all pair masks for each class.
       Class 0: 1080 pairs
       Class 1-3: up to 1080*1079 pairs total */
    int max_pairs = L * L;
    mask_t *pair_masks[4];
    int pair_counts[4] = {0, 0, 0, 0};

    for (int c = 0; c < 4; c++) {
        pair_masks[c] = (mask_t*)malloc(max_pairs * sizeof(mask_t));
    }

    for (int i = 0; i < L; i++) {
        mask_t mi = line_masks[i];
        int di = line_dirs[i];
        for (int j = 0; j < L; j++) {
            if (i == j) {
                pair_masks[0][pair_counts[0]++] = mi;
                continue;
            }
            mask_t mj = line_masks[j];
            if (mi & mj) {
                pair_masks[1][pair_counts[1]++] = mi | mj;
            } else {
                if (di == line_dirs[j]) {
                    pair_masks[2][pair_counts[2]++] = mi | mj;
                } else {
                    pair_masks[3][pair_counts[3]++] = mi | mj;
                }
            }
        }
    }

    int sizes[4] = {3, 5, 6, 6};
    mask_t reps[4];
    for (int c = 0; c < 4; c++) reps[c] = pair_masks[c][0];

    for (int a = 0; a < 4; a++) {
        mask_t repA = reps[a];
        long long countA = pair_counts[a];
        int sizeA = sizes[a];

        for (int b = 0; b < 4; b++) {
            int sizeB = sizes[b];
            long long dist[7] = {0, 0, 0, 0, 0, 0, 0};
            for (int m_idx = 0; m_idx < pair_counts[b]; m_idx++) {
                int r = popcount128(repA & pair_masks[b][m_idx]);
                dist[r]++;
            }
            for (int r = 0; r <= 6; r++) {
                if (dist[r]) {
                    int k = sizeA + sizeB - r;
                    A_k[k] += (i128)countA * dist[r];
                }
            }
        }
    }

    for (int c = 0; c < 4; c++) free(pair_masks[c]);
}

static i128 F_from_Ak(i128 *A_k, int n) {
    i128 tot = 0;
    for (int k = 0; k <= n; k++) {
        if (A_k[k]) {
            tot += A_k[k] * nCk(81 - k, n - k);
        }
    }
    return tot;
}

long long p818_native(void) {
    build_geometry();

    i128 A_k[13];
    memset(A_k, 0, sizeof(A_k));
    compute_Ak(A_k);

    i128 result = F_from_Ak(A_k, 12);
    return (long long)result;
}
