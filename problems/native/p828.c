// Project Euler 828: Numbers Challenge
// For each of 200 problems, find the minimum score (sum of used numbers)
// to reach the target using +, -, *, / with positive integer intermediates.
// Answer: sum(3^n * s_n) mod 1005075251.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;
typedef uint64_t u64;

#define MOD 1005075251LL

typedef struct {
    int target;
    int nums[6];
    int count;
} Problem;

static Problem problems[200];
static int num_problems;

static void read_data(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "Cannot open %s\n", path); exit(1); }
    char line[256];
    num_problems = 0;
    while (fgets(line, sizeof(line), f)) {
        char *p = line;
        int target = 0;
        while (*p >= '0' && *p <= '9') { target = target * 10 + (*p - '0'); p++; }
        problems[num_problems].target = target;
        problems[num_problems].count = 0;
        p++;  // skip ':'
        while (*p) {
            if (*p == '\n' || *p == '\r' || *p == '\0') break;
            int val = 0;
            while (*p >= '0' && *p <= '9') { val = val * 10 + (*p - '0'); p++; }
            if (*p == ',') p++;
            if (val > 0) problems[num_problems].nums[problems[num_problems].count++] = val;
        }
        num_problems++;
    }
    fclose(f);
}

// Dynamic array for values
typedef struct {
    int *vals;
    int count;
    int cap;
} DArr;

static void da_init(DArr *a) { a->cap = 16; a->vals = (int *)malloc(a->cap * sizeof(int)); a->count = 0; }
static void da_free(DArr *a) { free(a->vals); a->vals = NULL; a->count = 0; }
static void da_add(DArr *a, int v) {
    if (a->count >= a->cap) { a->cap *= 2; a->vals = (int *)realloc(a->vals, a->cap * sizeof(int)); }
    a->vals[a->count++] = v;
}

static int cmp_int(const void *a, const void *b) {
    int x = *(const int *)a, y = *(const int *)b;
    return (x > y) - (x < y);
}

// Dedup a sorted array in place, return new count.
static int dedup_sorted(int *arr, int n) {
    if (n <= 1) return n;
    int w = 1;
    for (int r = 1; r < n; r++) {
        if (arr[r] != arr[r - 1]) arr[w++] = arr[r];
    }
    return w;
}

// Cache: maps sorted-tuple key -> list of achievable values
typedef struct {
    u64 key;
    int *vals;
    int count;
    int used;
} CacheSlot;

#define CACHE_CAP (1 << 20)
static CacheSlot cache_slots[CACHE_CAP];

static void cache_reset(void) {
    memset(cache_slots, 0, sizeof(cache_slots));
}

static inline unsigned cache_hash(u64 key) {
    return (unsigned)(key * 1099511628211ULL);
}

static CacheSlot *cache_find(u64 key) {
    unsigned h = cache_hash(key) & (CACHE_CAP - 1);
    while (cache_slots[h].used) {
        if (cache_slots[h].key == key) return &cache_slots[h];
        h = (h + 1) & (CACHE_CAP - 1);
    }
    return NULL;
}

static CacheSlot *cache_insert(u64 key) {
    unsigned h = cache_hash(key) & (CACHE_CAP - 1);
    while (cache_slots[h].used) {
        if (cache_slots[h].key == key) return &cache_slots[h];
        h = (h + 1) & (CACHE_CAP - 1);
    }
    cache_slots[h].used = 1;
    cache_slots[h].key = key;
    cache_slots[h].vals = NULL;
    cache_slots[h].count = 0;
    return &cache_slots[h];
}

static u64 encode_key(int *nums, int len) {
    int sorted[6];
    memcpy(sorted, nums, len * sizeof(int));
    for (int i = 1; i < len; i++) {
        int key = sorted[i];
        int j = i - 1;
        while (j >= 0 && sorted[j] > key) { sorted[j+1] = sorted[j]; j--; }
        sorted[j+1] = key;
    }
    u64 k = 0;
    for (int i = 0; i < len; i++) k = k * 201 + (u64)(sorted[i] + 1);
    return k;
}

// Generate all achievable values from nums[0..len-1].
// Collects values (with duplicates) then sorts and dedups before caching.
static void generate_values(int *nums, int len, DArr *out) {
    u64 key = encode_key(nums, len);
    CacheSlot *slot = cache_find(key);
    if (slot) {
        for (int i = 0; i < slot->count; i++) da_add(out, slot->vals[i]);
        return;
    }

    DArr result;
    da_init(&result);

    if (len == 1) {
        da_add(&result, nums[0]);
    } else {
        int total = 1 << len;
        for (int mask = 1; mask < total - 1; mask++) {
            int complement = (~mask) & (total - 1);
            if (mask > complement) continue;

            int left[6], right[6];
            int nl = 0, nr = 0;
            for (int i = 0; i < len; i++) {
                if (mask & (1 << i)) left[nl++] = nums[i];
                else right[nr++] = nums[i];
            }

            DArr vl, vr;
            da_init(&vl); da_init(&vr);
            generate_values(left, nl, &vl);
            generate_values(right, nr, &vr);

            for (int i = 0; i < vl.count; i++) {
                for (int j = 0; j < vr.count; j++) {
                    int a = vl.vals[i], b = vr.vals[j];
                    da_add(&result, a + b);
                    da_add(&result, a * b);
                    if (a > b) da_add(&result, a - b);
                    if (b > a) da_add(&result, b - a);
                    if (b > 0 && a % b == 0) da_add(&result, a / b);
                    if (a > 0 && b % a == 0) da_add(&result, b / a);
                }
            }
            da_free(&vl); da_free(&vr);
        }

        // Sort and dedup
        qsort(result.vals, result.count, sizeof(int), cmp_int);
        result.count = dedup_sorted(result.vals, result.count);
    }

    // Store in cache
    CacheSlot *s = cache_insert(key);
    s->vals = (int *)malloc(result.count * sizeof(int));
    memcpy(s->vals, result.vals, result.count * sizeof(int));
    s->count = result.count;

    for (int i = 0; i < result.count; i++) da_add(out, result.vals[i]);
    da_free(&result);
}

static int min_score(int target, int *nums, int count) {
    int best = 0;
    int total = 1 << count;
    for (int mask = 1; mask < total; mask++) {
        int subset[6];
        int ns = 0;
        int sum = 0;
        for (int i = 0; i < count; i++) {
            if (mask & (1 << i)) {
                subset[ns++] = nums[i];
                sum += nums[i];
            }
        }
        if (best > 0 && sum >= best) continue;

        DArr vs;
        da_init(&vs);
        generate_values(subset, ns, &vs);

        // Binary search for target since values are sorted
        int lo = 0, hi = vs.count - 1, found = 0;
        while (lo <= hi) {
            int mid = (lo + hi) / 2;
            if (vs.vals[mid] == target) { found = 1; break; }
            if (vs.vals[mid] < target) lo = mid + 1;
            else hi = mid - 1;
        }
        if (found) {
            if (best == 0 || sum < best) best = sum;
        }
        da_free(&vs);
    }
    return best;
}

long long p828_native(void) {
    read_data("data/0828_number_challenges.txt");
    cache_reset();

    i64 total = 0;
    i64 pow3 = 3;
    for (int n = 0; n < num_problems; n++) {
        int s = min_score(problems[n].target, problems[n].nums, problems[n].count);
        total = (total + pow3 * s) % MOD;
        pow3 = (pow3 * 3) % MOD;
    }
    return total;
}
