/* Project Euler 949
 *
 * Compute G(20, 7, 1001001011).
 *
 * Algorithm:
 * 1. Build a binary tree DP computing "u" values for all 2^n binary strings
 * 2. Create histograms of u values (all and cold-only)
 * 3. Compute histogram^a and histogram^b via convolution
 * 4. Count pairs with sum < 0 and cold pairs with sum = 0
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;

#define MOD 1001001011LL

/* ---- ceil_div_pow2 ---- */

static i64 ceil_div_pow2(i64 x, int s) {
    if (s == 0) return x;
    i64 d = (i64)1 << s;
    if (x >= 0) return (x + d - 1) >> s;
    return -((-x) >> s);
}

/* ---- simplest_between ---- */

static i64 simplest_between(i64 u, i64 d, int e) {
    for (int m = 0; m <= e; m++) {
        int s = e - m;
        i64 p_min = (u >> s) + 1;
        i64 p_max = ceil_div_pow2(d, s) - 1;
        if (p_min <= p_max) {
            i64 p;
            if (p_min > 0) p = p_min;
            else if (p_max < 0) p = p_max;
            else p = 0;
            if (m > 0 && p != 0 && (p & 1) == 0) {
                if (p + 1 <= p_max && ((p + 1) & 1)) p += 1;
                else if (p - 1 >= p_min && ((p - 1) & 1)) p -= 1;
            }
            return p << s;
        }
    }
    return 0;
}

/* ---- compute_u_hot ---- */

static i64 *u_full_arr;
static int *hot_arr;

static void compute_u_hot(int n) {
    int e = n;
    i64 scale = (i64)1 << e;
    i64 total = ((i64)1 << (n + 1)) - 1;

    i64 *dp_u = (i64 *)calloc(total, sizeof(i64));
    i64 *dp_d = (i64 *)calloc(total, sizeof(i64));

    /* Root's children at index 1, 2 */
    dp_u[1] = scale;  dp_d[1] = scale;    /* left: 1/1 */
    dp_u[2] = -scale; dp_d[2] = -scale;   /* right: -1/-1 */

    u_full_arr = (i64 *)calloc((i64)1 << n, sizeof(i64));
    hot_arr = (int *)calloc((i64)1 << n, sizeof(int));

    for (int length = 2; length <= n; length++) {
        i64 size = (i64)1 << length;
        i64 start = ((i64)1 << length) - 1;

        for (i64 bits = 0; bits < size; bits++) {
            /* u_raw = max over suffixes of dp_d */
            i64 u_raw = INT64_MIN;
            for (int s_len = 1; s_len < length; s_len++) {
                i64 suf = bits & (((i64)1 << s_len) - 1);
                i64 cand = dp_d[((i64)1 << s_len) - 1 + suf];
                if (cand > u_raw) u_raw = cand;
            }

            /* d_raw = min over prefixes of dp_u */
            i64 d_raw = INT64_MAX;
            for (int p_len = 1; p_len < length; p_len++) {
                i64 pre = bits >> (length - p_len);
                i64 cand = dp_u[((i64)1 << p_len) - 1 + pre];
                if (cand < d_raw) d_raw = cand;
            }

            i64 idx = start + bits;
            if (u_raw < d_raw) {
                i64 x = simplest_between(u_raw, d_raw, e);
                dp_u[idx] = x;
                dp_d[idx] = x;
                if (length == n) {
                    u_full_arr[bits] = x;
                    hot_arr[bits] = 0;
                }
            } else {
                dp_u[idx] = u_raw;
                dp_d[idx] = d_raw;
                if (length == n) {
                    u_full_arr[bits] = u_raw;
                    hot_arr[bits] = 1;
                }
            }
        }
    }

    free(dp_u);
    free(dp_d);
}

/* ---- Histogram (sorted array of (val, cnt) pairs) ---- */

typedef struct { i64 val; i64 cnt; } HistEntry;

typedef struct {
    HistEntry *entries;
    int size;
    int cap;
} Histogram;

static void hist_init(Histogram *h) {
    h->cap = 16;
    h->size = 0;
    h->entries = (HistEntry *)malloc(h->cap * sizeof(HistEntry));
}

static void hist_free(Histogram *h) {
    free(h->entries);
    h->entries = NULL;
    h->size = 0;
    h->cap = 0;
}

static int hist_find(Histogram *h, i64 val) {
    int lo = 0, hi = h->size;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (h->entries[mid].val < val) lo = mid + 1;
        else hi = mid;
    }
    if (lo < h->size && h->entries[lo].val == val) return lo;
    return -1;
}

static void hist_insert_sorted(Histogram *h, i64 val, i64 cnt) {
    if (cnt == 0) return;
    int pos = hist_find(h, val);
    if (pos >= 0) {
        h->entries[pos].cnt = (h->entries[pos].cnt + cnt) % MOD;
        if (h->entries[pos].cnt == 0) {
            /* Remove zero entry */
            memmove(&h->entries[pos], &h->entries[pos + 1],
                    (h->size - pos - 1) * sizeof(HistEntry));
            h->size--;
        }
        return;
    }
    /* Insert at sorted position */
    int lo = 0, hi = h->size;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (h->entries[mid].val < val) lo = mid + 1;
        else hi = mid;
    }
    if (h->size >= h->cap) {
        h->cap *= 2;
        h->entries = (HistEntry *)realloc(h->entries, h->cap * sizeof(HistEntry));
    }
    memmove(&h->entries[lo + 1], &h->entries[lo],
            (h->size - lo) * sizeof(HistEntry));
    h->entries[lo].val = val;
    h->entries[lo].cnt = cnt % MOD;
    if (h->entries[lo].cnt == 0) return; /* shouldn't happen since cnt != 0 */
    h->size++;
}

static Histogram hist_from_values(i64 *values, int n) {
    Histogram h;
    hist_init(&h);
    for (int i = 0; i < n; i++) {
        int pos = hist_find(&h, values[i]);
        if (pos >= 0) {
            h.entries[pos].cnt = (h.entries[pos].cnt + 1) % MOD;
            if (h.entries[pos].cnt == 0) {
                memmove(&h.entries[pos], &h.entries[pos + 1],
                        (h.size - pos - 1) * sizeof(HistEntry));
                h.size--;
            }
        } else {
            hist_insert_sorted(&h, values[i], 1);
        }
    }
    return h;
}

static Histogram hist_convolve(Histogram *a, Histogram *b) {
    Histogram out;
    hist_init(&out);

    /* Iterate over smaller */
    Histogram *smaller = a, *larger = b;
    if (a->size > b->size) { smaller = b; larger = a; }

    for (int i = 0; i < smaller->size; i++) {
        for (int j = 0; j < larger->size; j++) {
            i64 key = smaller->entries[i].val + larger->entries[j].val;
            i64 cnt = (i64)((__int128)smaller->entries[i].cnt * larger->entries[j].cnt % MOD);
            hist_insert_sorted(&out, key, cnt);
        }
    }
    return out;
}

static Histogram hist_pow(Histogram *h, int t) {
    if (t == 0) {
        Histogram r;
        hist_init(&r);
        hist_insert_sorted(&r, 0, 1);
        return r;
    }
    Histogram result;
    hist_init(&result);
    /* Copy h into result */
    for (int i = 0; i < h->size; i++)
        hist_insert_sorted(&result, h->entries[i].val, h->entries[i].cnt);

    for (int i = 1; i < t; i++) {
        Histogram next = hist_convolve(&result, h);
        hist_free(&result);
        result = next;
    }
    return result;
}

/* ---- Counting functions ---- */

static i64 count_sum_lt_zero(Histogram *a, Histogram *b) {
    /* Sort b by val (already sorted), compute prefix sums of counts */
    i64 *b_sums = (i64 *)malloc(b->size * sizeof(i64));
    i64 *pref = (i64 *)malloc((b->size + 1) * sizeof(i64));

    for (int i = 0; i < b->size; i++) b_sums[i] = b->entries[i].val;
    pref[0] = 0;
    for (int i = 0; i < b->size; i++)
        pref[i + 1] = (pref[i] + b->entries[i].cnt) % MOD;

    i64 ans = 0;
    for (int i = 0; i < a->size; i++) {
        i64 sa = a->entries[i].val;
        i64 ca = a->entries[i].cnt;
        /* bisect_left for -sa in b_sums */
        int lo = 0, hi = b->size;
        while (lo < hi) {
            int mid = lo + (hi - lo) / 2;
            if (b_sums[mid] < -sa) lo = mid + 1;
            else hi = mid;
        }
        ans = (ans + (i64)((__int128)ca * pref[lo] % MOD)) % MOD;
    }

    free(b_sums);
    free(pref);
    return ans;
}

static i64 count_sum_eq_zero(Histogram *a, Histogram *b) {
    Histogram *smaller = a, *larger = b;
    if (a->size > b->size) { smaller = b; larger = a; }

    i64 ans = 0;
    for (int i = 0; i < smaller->size; i++) {
        i64 s = smaller->entries[i].val;
        i64 ca = smaller->entries[i].cnt;
        int pos = hist_find(larger, -s);
        if (pos >= 0) {
            ans = (ans + (i64)((__int128)ca * larger->entries[pos].cnt % MOD)) % MOD;
        }
    }
    return ans;
}

/* ---- Main solver ---- */

long long p949_native(void) {
    int n = 20;
    int k = 7;

    compute_u_hot(n);

    int total = (i64)1 << n;

    /* u_hist from all u_full values */
    Histogram u_hist = hist_from_values(u_full_arr, total);

    /* cold_hist from cold values only */
    i64 *cold_values = (i64 *)malloc(total * sizeof(i64));
    int n_cold = 0;
    for (int i = 0; i < total; i++) {
        if (hot_arr[i] == 0) cold_values[n_cold++] = u_full_arr[i];
    }
    Histogram cold_hist = hist_from_values(cold_values, n_cold);

    int a_pow = k / 2;  /* 3 */
    int b_pow = k - a_pow;  /* 4 */

    Histogram dist_a = hist_pow(&u_hist, a_pow);
    Histogram dist_b = hist_pow(&u_hist, b_pow);

    i64 neg = count_sum_lt_zero(&dist_a, &dist_b);

    Histogram cold_a = hist_pow(&cold_hist, a_pow);
    Histogram cold_b = hist_pow(&cold_hist, b_pow);

    i64 zero_cold = count_sum_eq_zero(&cold_a, &cold_b);

    i64 result = (neg + zero_cold) % MOD;
    if (result < 0) result += MOD;

    /* Cleanup */
    hist_free(&u_hist);
    hist_free(&cold_hist);
    hist_free(&dist_a);
    hist_free(&dist_b);
    hist_free(&cold_a);
    hist_free(&cold_b);
    free(cold_values);
    free(u_full_arr);
    free(hot_arr);

    return result;
}
