// Project Euler 846: Magic Bracelets
// F(10^6) using Gaussian integer vectors and path merging.
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef long long i64;
typedef unsigned long long u64;

#define LIMIT 1000000

// Vector storage
static struct { i64 x, y; char valid; } vec[LIMIT + 1];

// Square lookup: sq_root[n] = sqrt(n) if n is a perfect square, else 0
static int sq_root[LIMIT + 1];

// Sieve
static char is_prime[LIMIT + 1];

// Extended GCD for modular inverse
static i64 ext_gcd(i64 a, i64 b, i64 *x, i64 *y) {
    if (b == 0) { *x = 1; *y = 0; return a; }
    i64 x1, y1;
    i64 g = ext_gcd(b, a % b, &x1, &y1);
    *x = y1;
    *y = x1 - (a / b) * y1;
    return g;
}

static i64 mod_inverse(i64 a, i64 m) {
    i64 x, y;
    i64 g = ext_gcd(((a % m) + m) % m, m, &x, &y);
    if (g != 1) return -1;
    return ((x % m) + m) % m;
}

static void canonical(i64 *x, i64 *y) {
    i64 a = *x < 0 ? -*x : *x;
    i64 b = *y < 0 ? -*y : *y;
    if (b > a) { i64 t = a; a = b; b = t; }
    *x = a; *y = b;
}

static void build_vectors(int limit) {
    vec[1].x = 1; vec[1].y = 0; vec[1].valid = 1;
    vec[2].x = 1; vec[2].y = 1; vec[2].valid = 1;

    // Build square root lookup
    memset(sq_root, 0, sizeof(sq_root));
    for (int i = 0; (i64)i * i <= limit; i++)
        sq_root[i * i] = i;

    // Sieve
    memset(is_prime, 1, sizeof(is_prime));
    is_prime[0] = is_prime[1] = 0;
    for (int p = 2; (i64)p * p <= limit; p++) {
        if (is_prime[p]) {
            for (int m = p * p; m <= limit; m += p)
                is_prime[m] = 0;
        }
    }

    // For each prime p = 1 (mod 4), find base representation
    for (int p = 5; p <= limit; p++) {
        if (!is_prime[p] || p % 4 != 1) continue;

        i64 bx = 0, by = 0;
        for (int x = 1; (i64)x * x <= p; x++) {
            int rem = p - x * x;
            if (rem >= 0 && rem <= limit && sq_root[rem]) {
                bx = x; by = sq_root[rem];
                canonical(&bx, &by);
                break;
            }
        }
        if (bx == 0) continue;

        i64 value = p;
        i64 x = bx, y = by;
        while (value <= limit) {
            i64 cx = x, cy = y;
            canonical(&cx, &cy);
            vec[value].x = cx; vec[value].y = cy; vec[value].valid = 1;

            i64 doubled = 2 * value;
            if (doubled <= limit) {
                i64 dx = cx - cy, dy = cx + cy;
                canonical(&dx, &dy);
                vec[doubled].x = dx; vec[doubled].y = dy; vec[doubled].valid = 1;
            }

            i64 nx = x * bx - y * by;
            i64 ny = x * by + y * bx;
            x = nx; y = ny;
            value *= p;
        }
    }
}

// Edge hash map
#define EDGE_HT_SIZE (1 << 21)
#define EDGE_HT_MASK (EDGE_HT_SIZE - 1)

typedef struct {
    i64 key;
    i64 count;
    i64 sum;
    char used;
} EdgeEntry;

static EdgeEntry *edge_ht;

static i64 encode_edge(i64 a, i64 b) {
    if (a < b) { i64 t = a; a = b; b = t; }
    return a * (i64)(LIMIT + 1) + b;
}

static unsigned hash_key(i64 key) {
    u64 h = (u64)key;
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 33;
    return (unsigned)(h & EDGE_HT_MASK);
}

static EdgeEntry *edge_find(i64 key) {
    unsigned h = hash_key(key);
    while (edge_ht[h].used) {
        if (edge_ht[h].key == key) return &edge_ht[h];
        h = (h + 1) & EDGE_HT_MASK;
    }
    return &edge_ht[h];
}

static EdgeEntry *edge_lookup(i64 key) {
    unsigned h = hash_key(key);
    while (edge_ht[h].used) {
        if (edge_ht[h].key == key) return &edge_ht[h];
        h = (h + 1) & EDGE_HT_MASK;
    }
    return NULL;
}

// Parent labels
static void parent_labels(i64 label, i64 x, i64 y, i64 *parents, int *np) {
    *np = 0;
    if (label == 2) {
        parents[(*np)++] = 1;
        return;
    }
    if (y == 0) return;

    i64 beta = mod_inverse(y, x);
    if (beta < 0) return;
    i64 alpha = (1 - beta * y) / x;

    i64 q1x = beta < 0 ? -beta : beta;
    i64 q1y = alpha < 0 ? -alpha : alpha;
    i64 q2x = x - q1x;
    i64 q2y = y - q1y;

    i64 n1 = q1x * q1x + q1y * q1y;
    i64 n2 = q2x * q2x + q2y * q2y;

    if (n1 <= LIMIT && vec[n1].valid) {
        parents[(*np)++] = n1;
    }
    if (n2 <= LIMIT && vec[n2].valid && n2 != n1) {
        parents[(*np)++] = n2;
    }
    // Sort descending
    if (*np == 2 && parents[0] < parents[1]) {
        i64 t = parents[0]; parents[0] = parents[1]; parents[1] = t;
    }
}

long long p846_native(void) {
    build_vectors(LIMIT);

    // Edge hash map
    edge_ht = calloc(EDGE_HT_SIZE, sizeof(EdgeEntry));
    if (!edge_ht) return -1;

    i64 total = 0;

    // Process labels in decreasing order
    for (int label = LIMIT; label >= 2; label--) {
        if (!vec[label].valid) continue;

        i64 x = vec[label].x, y = vec[label].y;
        i64 parents[2];
        int np;
        parent_labels(label, x, y, parents, &np);

        if (np == 0) continue;

        for (int pi = 0; pi < np; pi++) {
            i64 parent = parents[pi];
            i64 key = encode_edge(label, parent);
            EdgeEntry *e = edge_find(key);
            i64 count = 0, isum = 0;
            if (e->used) {
                count = e->count;
                isum = e->sum;
            }

            total += count * (label + parent) + isum;

            if (!e->used) {
                e->used = 1;
                e->key = key;
            }
            e->count = count + 1;
            e->sum = isum;
        }

        if (np == 2) {
            i64 a = parents[0], b = parents[1];
            i64 key_a = encode_edge(label, a);
            i64 key_b = encode_edge(label, b);
            EdgeEntry *ea = edge_lookup(key_a);
            EdgeEntry *eb = edge_lookup(key_b);

            // ea and eb must exist since we just inserted them above
            if (!ea || !eb) continue;

            i64 count_a = ea->count, sum_a = ea->sum;
            i64 count_b = eb->count, sum_b = eb->sum;

            i64 merged_count = count_a * count_b;
            i64 merged_sum = merged_count * label + count_a * sum_b + count_b * sum_a;

            i64 lower_key = encode_edge(a, b);
            EdgeEntry *le = edge_find(lower_key);
            i64 old_count = 0, old_sum = 0;
            if (le->used) {
                old_count = le->count;
                old_sum = le->sum;
            }
            if (!le->used) {
                le->used = 1;
                le->key = lower_key;
            }
            le->count = old_count + merged_count;
            le->sum = old_sum + merged_sum;
        }
    }

    free(edge_ht);
    return total;
}
