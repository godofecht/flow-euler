// Project Euler 269: inclusion-exclusion over integer roots -1 through -9.
// A state contains up to nine signed bytes, so use an unsigned __int128 key.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

enum { CAP = 200003 };

typedef struct {
    unsigned __int128 *keys;
    int64_t *vals;
    unsigned char *used;
} Map;

static unsigned long hash_key(unsigned __int128 key) {
    uint64_t lo = (uint64_t)key, hi = (uint64_t)(key >> 64);
    uint64_t h = lo ^ (hi + 0x9e3779b97f4a7c15ULL + (lo << 6) + (lo >> 2));
    return (unsigned long)(h % CAP);
}

static long find_slot(const Map *map, unsigned __int128 key) {
    long h = (long)hash_key(key);
    while (map->used[h] && map->keys[h] != key) if (++h == CAP) h = 0;
    return h;
}

static void add(Map *map, unsigned __int128 key, int64_t value) {
    long h = find_slot(map, key);
    if (!map->used[h]) {
        map->used[h] = 1;
        map->keys[h] = key;
        map->vals[h] = value;
    } else map->vals[h] += value;
}

static void clear(Map *map) {
    for (long i = 0; i < CAP; ++i) map->used[i] = 0;
}

static unsigned __int128 encode(const int64_t *coeffs, int m) {
    unsigned __int128 key = 0;
    for (int i = 0; i < m; ++i) key = (key << 8) | (unsigned char)(coeffs[i] + 128);
    return key;
}

static void decode(unsigned __int128 key, int64_t *coeffs, int m) {
    for (int i = m - 1; i >= 0; --i) {
        coeffs[i] = (int64_t)(key & 255) - 128;
        key >>= 8;
    }
}

static int64_t count_for(const int64_t *roots, int m, int length) {
    Map a = { calloc(CAP, sizeof(*a.keys)), calloc(CAP, sizeof(*a.vals)), calloc(CAP, 1) };
    Map b = { calloc(CAP, sizeof(*b.keys)), calloc(CAP, sizeof(*b.vals)), calloc(CAP, 1) };
    int64_t coeffs[9] = {0}, next[9];
    unsigned __int128 zero = encode(coeffs, m);
    add(&a, zero, 1);
    for (int pos = 0; pos < length; ++pos) {
        clear(&b);
        int lo = (pos == 0 || pos == length - 1) ? 1 : 0;
        for (long h = 0; h < CAP; ++h) if (a.used[h]) {
            decode(a.keys[h], coeffs, m);
            for (int d = lo; d <= 9; ++d) {
                int ok = 1;
                for (int j = 0; j < m; ++j) {
                    int64_t sum = coeffs[j] + d;
                    if (sum % roots[j]) { ok = 0; break; }
                    next[j] = -sum / roots[j];
                }
                if (ok) add(&b, encode(next, m), a.vals[h]);
            }
        }
        Map tmp = a; a = b; b = tmp;
    }
    long h = find_slot(&a, zero);
    int64_t answer = a.used[h] ? a.vals[h] : 0;
    free(a.keys); free(a.vals); free(a.used);
    free(b.keys); free(b.vals); free(b.used);
    return answer;
}

int64_t p269_native(void) {
    const int64_t base = 1000000000000000LL;
    int64_t roots[9], total = 0;
    for (int mask = 1; mask < (1 << 9); ++mask) {
        int m = 0;
        for (int i = 0; i < 9; ++i) if (mask & (1 << i)) roots[m++] = i + 1;
        int64_t subtotal = 0;
        for (int len = 1; len <= 16; ++len) subtotal += count_for(roots, m, len);
        total += (m & 1) ? subtotal : -subtotal;
    }
    return base + total;
}
