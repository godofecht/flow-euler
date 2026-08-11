/* Project Euler 963: A Gated Town
 *
 * F(N) for N = 10^5.
 *
 * For each n in 1..N, compute a key (L_tuple, parity_of_2s, B_value).
 * Count occurrences of each key, then combine pairs and sum squares / 4.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

static double *B;

static void build_B(int m) {
    B = calloc(m + 1, sizeof(double));
    for (int i = 1; i <= m; i++) {
        double h = B[i >> 1];
        if (i & 1) B[i] = h + 1;
        else if (h < 2) B[i] = h / 2;
        else B[i] = h - 1;
    }
}

/* Key: (L_tuple, parity, B_val)
 * Encode as a sortable byte string:
 *   1 byte: len
 *   len * 4 bytes: L values (int32)
 *   4 bytes: parity (int32)
 *   8 bytes: bval (double)
 * Total: up to 1 + 16*4 + 4 + 8 = 77 bytes
 */

typedef struct {
    int len;
    int L[24];
    int parity;
    double bval;
} Key;

/* Compare keys for sorting (lexicographic) */
static int cmp_key(const void *a, const void *b) {
    const Key *ka = (const Key *)a;
    const Key *kb = (const Key *)b;
    if (ka->len != kb->len) return ka->len < kb->len ? -1 : 1;
    for (int i = 0; i < ka->len; i++) {
        if (ka->L[i] != kb->L[i]) return ka->L[i] < kb->L[i] ? -1 : 1;
    }
    if (ka->parity != kb->parity) return ka->parity < kb->parity ? -1 : 1;
    if (ka->bval < kb->bval) return -1;
    if (ka->bval > kb->bval) return 1;
    return 0;
}

static int keys_equal(const Key *a, const Key *b) {
    if (a->len != b->len) return 0;
    for (int i = 0; i < a->len; i++)
        if (a->L[i] != b->L[i]) return 0;
    if (a->parity != b->parity) return 0;
    return a->bval == b->bval;
}

static void compute_key(i64 n, Key *out) {
    char s[32];
    int slen = 0;
    i64 x = n;
    while (x) {
        s[slen++] = (char)(x % 3) + '0';
        x /= 3;
    }
    for (int i = 0; i < slen / 2; i++) {
        char t = s[i]; s[i] = s[slen-1-i]; s[slen-1-i] = t;
    }
    if (slen == 0) { s[0] = '0'; slen = 1; }

    int l = -1;
    for (int j = 0; j < slen; j++) {
        if (s[j] == '1') break;
        if (s[j] < '1') l = j;
    }

    out->len = 0;
    if (l > -1) {
        int k = 0;
        for (int j = l - 1; j >= 0; j--) {
            if (s[j] < '1') k++;
            else {
                out->L[out->len++] = k;
            }
        }
    }

    int cnt2 = 0;
    for (int j = 0; j < slen; j++)
        if (s[j] == '2') cnt2++;
    out->parity = cnt2 & 1;

    char u[32];
    int ulen = 0;
    for (int j = 0; j < slen; j++)
        if (s[j] < '2') u[ulen++] = s[j];
    if (ulen == 0) {
        out->bval = 0.0;
    } else {
        i64 val = 0;
        for (int j = 0; j < ulen; j++)
            val = val * 2 + (u[j] - '0');
        out->bval = B[val];
    }
    for (int i = out->len; i < 24; i++) out->L[i] = 0;
}

static Key combine_key(const Key *a, const Key *b) {
    Key out;
    out.len = a->len + b->len;
    int idx = 0;
    for (int i = 0; i < a->len; i++) out.L[idx++] = a->L[i];
    for (int i = 0; i < b->len; i++) out.L[idx++] = b->L[i];
    /* Sort ascending */
    for (int i = 0; i < out.len; i++)
        for (int j = i + 1; j < out.len; j++)
            if (out.L[i] > out.L[j]) {
                int t = out.L[i]; out.L[i] = out.L[j]; out.L[j] = t;
            }
    out.parity = a->parity ^ b->parity;
    out.bval = a->bval + b->bval;
    for (int i = out.len; i < 24; i++) out.L[i] = 0;
    return out;
}

long long p963_native(void) {
    i64 N = 100000;
    int m = 1 << (int)ceil(log((double)(N > 1 ? N : 1)) / log(3.0));
    if (m < 1) m = 1;
    build_B(m);

    /* Build all keys for n=1..N */
    Key *all_keys = malloc(N * sizeof(Key));
    for (i64 n = 1; n <= N; n++) {
        compute_key(n, &all_keys[n - 1]);
    }

    /* Sort and count unique keys */
    qsort(all_keys, N, sizeof(Key), cmp_key);

    /* Extract unique keys and their counts */
    i64 max_unique = N;
    Key *unique_keys = malloc(max_unique * sizeof(Key));
    i64 *unique_counts = malloc(max_unique * sizeof(i64));
    i64 nkeys = 0;

    i64 i = 0;
    while (i < N) {
        i64 j = i + 1;
        while (j < N && keys_equal(&all_keys[i], &all_keys[j])) j++;
        unique_keys[nkeys] = all_keys[i];
        unique_counts[nkeys] = j - i;
        nkeys++;
        i = j;
    }

    free(all_keys);

    /* Generate all combined keys with their counts */
    i64 total_combos = nkeys * nkeys;
    /* Each combo produces a key and a count. Sort and merge. */
    /* Use a struct to hold key + count for sorting */
    typedef struct {
        Key key;
        i64 count;
    } Combo;
    Combo *combos = malloc(total_combos * sizeof(Combo));
    i64 ci = 0;
    for (i64 ai = 0; ai < nkeys; ai++) {
        for (i64 bi = 0; bi < nkeys; bi++) {
            combos[ci].key = combine_key(&unique_keys[ai], &unique_keys[bi]);
            combos[ci].count = unique_counts[ai] * (unique_counts[bi] + (ai == bi ? 1 : 0));
            ci++;
        }
    }

    free(unique_keys);
    free(unique_counts);

    /* Sort combos by key */
    qsort(combos, total_combos, sizeof(Combo), cmp_key);

    /* Merge identical keys and sum counts */
    i128 ans = 0;
    i64 k = 0;
    while (k < total_combos) {
        i64 j = k + 1;
        while (j < total_combos && keys_equal(&combos[k].key, &combos[j].key)) j++;
        i128 v = 0;
        for (i64 x = k; x < j; x++) v += combos[x].count;
        ans += v * v;
        k = j;
    }
    ans /= 4;

    free(combos);
    free(B);
    return (long long)ans;
}
