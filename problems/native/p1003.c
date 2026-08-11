// Project Euler 1003: Lonely Singles - S(80).
// Meet-in-the-middle enumeration with gap-3 constraint.
// Adapted from the Python solver by cirosantilli.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    long long boundary;
    long long tail;
    long long value;
    int prefix_mask;
} Entry;

static int cmp_entry(const void *a, const void *b) {
    const Entry *ea = (const Entry *)a;
    const Entry *eb = (const Entry *)b;
    if (ea->boundary != eb->boundary) return ea->boundary < eb->boundary ? -1 : 1;
    if (ea->prefix_mask != eb->prefix_mask) return ea->prefix_mask < eb->prefix_mask ? -1 : 1;
    if (ea->tail != eb->tail) return ea->tail < eb->tail ? -1 : 1;
    return 0;
}

// Coefficients
static long long h[84];
static long long boundary[80];
static long long tail[80];
static long long value[80];

static Entry *right_entries;
static long long *right_prefix; // prefix sums of value
static int right_count;

static int limit;
static int split;

// Generate right half configurations
static void gen_right(int pos, int prev1, int prev2,
                      long long bs, long long ts, long long vs, int pm) {
    if (pos == limit) {
        right_entries[right_count].boundary = bs;
        right_entries[right_count].tail = ts;
        right_entries[right_count].value = vs;
        right_entries[right_count].prefix_mask = pm;
        right_count++;
        return;
    }
    // Don't place at pos
    gen_right(pos + 1, 0, prev1, bs, ts, vs, pm);
    if (prev1 || prev2) return;
    // Place at pos
    int npm = pm;
    if (pos == split) npm |= 1;
    else if (pos == split + 1) npm |= 2;
    gen_right(pos + 1, 1, prev1, bs + boundary[pos], ts + tail[pos], vs + value[pos], npm);
}

static int allowed_masks[4][4];
static int allowed_count[4];

static void init_allowed(void) {
    for (int mask = 0; mask < 4; mask++) {
        allowed_count[mask] = 0;
        for (int prefix = 0; prefix < 4; prefix++) {
            if ((mask & 1) && (prefix & 3)) continue;
            if ((mask & 2) && (prefix & 1)) continue;
            allowed_masks[mask][allowed_count[mask]++] = prefix;
        }
    }
}

static long long total_sum;

static void gen_left(int pos, int prev1, int prev2,
                     long long bs, long long ts, long long vs) {
    if (pos == split) {
        int last_mask = prev1 | (prev2 << 1);
        long long needed_boundary = -bs;
        long long needed_tail = -ts;
        for (int i = 0; i < allowed_count[last_mask]; i++) {
            int pm = allowed_masks[last_mask][i];
            // Binary search for entries with (needed_boundary, pm) and tail >= needed_tail
            Entry key;
            key.boundary = needed_boundary;
            key.prefix_mask = pm;
            key.tail = needed_tail;
            // Find first entry >= key
            int lo = 0, hi = right_count;
            while (lo < hi) {
                int mid = (lo + hi) / 2;
                if (cmp_entry(&right_entries[mid], &key) < 0)
                    lo = mid + 1;
                else
                    hi = mid;
            }
            int start = lo;
            // Find first entry with (needed_boundary, pm) and tail > needed_tail
            // Actually we need tail >= needed_tail, which is what we found (start)
            // Find end of group (needed_boundary, pm): first entry with
            // (boundary, prefix_mask) > (needed_boundary, pm)
            key.boundary = needed_boundary;
            key.prefix_mask = pm + 1;
            key.tail = -9223372036854775807LL;
            lo = start; hi = right_count;
            while (lo < hi) {
                int mid = (lo + hi) / 2;
                if (cmp_entry(&right_entries[mid], &key) < 0)
                    lo = mid + 1;
                else
                    hi = mid;
            }
            int end = lo;
            int count = end - start;
            if (count > 0) {
                total_sum += (long long)count * vs + right_prefix[end] - right_prefix[start];
            }
        }
        return;
    }
    // Don't place at pos
    gen_left(pos + 1, 0, prev1, bs, ts, vs);
    if (prev1 || prev2) return;
    // Place at pos
    long long add_b = 0, add_t = 0, add_v = 0;
    if (pos == 0) { add_v = 1; }
    else if (pos == 1) { add_b = -1; }
    else if (pos == 2) { add_b = 1; add_t = -1; add_v = -2; }
    else { add_b = boundary[pos]; add_t = tail[pos]; add_v = value[pos]; }
    gen_left(pos + 1, 1, prev1, bs + add_b, ts + add_t, vs + add_v);
}

long long p1003_native(void) {
    limit = 80;
    split = limit / 2;

    // Compute coefficients
    h[1] = 0; h[2] = 0; h[3] = 1;
    for (int i = 4; i <= limit; i++)
        h[i] = 2 * h[i - 3] - h[i - 2];
    for (int pos = 3; pos < limit; pos++) {
        boundary[pos] = h[pos] - 3 * h[pos - 1] + 2 * h[pos - 2];
        tail[pos] = h[pos - 1] - 2 * h[pos - 2];
        value[pos] = 2 * (tail[pos] + h[pos]);
    }

    init_allowed();

    // Allocate right entries (generous upper bound)
    right_entries = (Entry *)malloc(8000000 * sizeof(Entry));
    right_prefix = (long long *)malloc(8000001 * sizeof(long long));
    right_count = 0;

    gen_right(split, 0, 0, 0, 0, 0, 0);

    // Sort right entries
    qsort(right_entries, right_count, sizeof(Entry), cmp_entry);

    // Build prefix sums of value
    right_prefix[0] = 0;
    for (int i = 0; i < right_count; i++)
        right_prefix[i + 1] = right_prefix[i] + right_entries[i].value;

    // Generate left half and accumulate
    total_sum = 0;
    gen_left(0, 0, 0, 0, 0, 0);

    free(right_entries);
    free(right_prefix);

    return total_sum;
}
