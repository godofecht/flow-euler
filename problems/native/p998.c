// Project Euler 998: Minimum Bounding Square of Triangles
// Pythagorean partners approach with exact integer arithmetic.
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

#define LIMIT 1000000

// Integer square root.
static i64 isqrt_i64(i64 n) {
    if (n <= 0) return 0;
    i64 r = (i64)sqrt((double)n);
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r;
}

static i64 gcd_i64(i64 a, i64 b) {
    while (b) { i64 t = a % b; a = b; b = t; }
    return a;
}

// Partners storage: for each m, a list of x values.
// We use a single allocation with per-m offsets.
static i64 *partner_data;
static int *partner_count;
static int *partner_cap;
static i64 **partner_ptr;

static void partners_add(int m, i64 x) {
    if (partner_count[m] >= partner_cap[m]) {
        int newcap = partner_cap[m] * 2 + 4;
        i64 *newptr = realloc(partner_ptr[m], newcap * sizeof(i64));
        partner_ptr[m] = newptr;
        partner_cap[m] = newcap;
    }
    partner_ptr[m][partner_count[m]++] = x;
}

static void build_partners(int limit) {
    partner_count = calloc(limit + 1, sizeof(int));
    partner_cap = calloc(limit + 1, sizeof(int));
    partner_ptr = calloc(limit + 1, sizeof(i64*));

    for (i64 r = 2; r * r <= 2 * (i64)limit; r++) {
        i64 rr = r * r;
        for (i64 s = 1; s < r; s++) {
            if (((r - s) & 1) == 0) continue;
            if (gcd_i64(r, s) != 1) continue;
            i64 a = rr - s * s;
            i64 b = 2 * r * s;
            i64 m = (a > b) ? a : b;
            i64 x = (a > b) ? b : a;
            if (m > limit) continue;
            for (i64 km = m; km <= limit; km += m) {
                partners_add((int)km, (km / m) * x);
            }
        }
    }

    // Sort each partner list.
    for (int m = 1; m <= limit; m++) {
        if (partner_count[m] > 1) {
            // Simple insertion sort (lists are small).
            for (int i = 1; i < partner_count[m]; i++) {
                i64 key = partner_ptr[m][i];
                int j = i - 1;
                while (j >= 0 && partner_ptr[m][j] > key) {
                    partner_ptr[m][j + 1] = partner_ptr[m][j];
                    j--;
                }
                partner_ptr[m][j + 1] = key;
            }
        }
    }
}

// Check that a triangle's minimum bounding square has the given side.
static int is_minimum_square(i64 sides[3], i64 twice_area, i64 square_side) {
    i64 m = square_side;
    i64 m2 = m * m;
    i64 d_area = twice_area;
    int has_equal = 0;

    // Candidate 1: square side parallel to a triangle side.
    for (int i = 0; i < 3; i++) {
        i64 d = sides[i];
        i64 e = sides[(i + 1) % 3];
        i64 f = sides[(i + 2) % 3];
        i64 den = 2 * d;
        i64 t_num = d * d + e * e - f * f;
        i64 d_num = d * den;
        i64 mx = d_num, mn = d_num;
        if (t_num > mx) mx = t_num;
        if (t_num < mn) mn = t_num;
        if (0 > mx) mx = 0;
        if (0 < mn) mn = 0;
        i64 width_num = mx - mn;

        i64 width_cmp = width_num - m * den;
        i64 height_cmp = d_area - m * d;
        if (width_cmp < 0 && height_cmp < 0)
            return 0;
        if (width_cmp <= 0 && height_cmp <= 0 && (width_cmp == 0 || height_cmp == 0))
            has_equal = 1;
    }

    // Candidate 2: all four sides of the square are touched.
    for (int i = 0; i < 3; i++) {
        i64 r = sides[i];
        i64 p = sides[(i + 1) % 3];
        i64 q = sides[(i + 2) % 3];
        i64 k_num = p * p + q * q - r * r;  // 2K
        if (k_num <= 0) continue;
        i64 r_den_part = p * p + q * q - 2 * d_area;
        if (r_den_part <= 0) continue;

        i128 num = (i128)k_num * k_num;  // 4B^2 numerator
        i128 den = (i128)4 * r_den_part;  // B^2 denominator

        i64 p2 = p * p;
        i64 q2 = q * q;

        if (num < (i128)d_area * den) continue;
        if (num > (i128)p2 * den || num > (i128)q2 * den) continue;
        if ((i128)p2 * den > 2 * num || (i128)q2 * den > 2 * num) continue;

        i128 target = (i128)m2 * den;
        if (num < target)
            return 0;
        if (num == target)
            has_equal = 1;
    }

    return has_equal;
}

// Hash set for sorted side triples.
// Each side fits in 21 bits (up to ~1.5M). Key = (a<<42)|(b<<21)|c.
#define HT_SIZE (1 << 24)
#define HT_MASK (HT_SIZE - 1)

static u64 *ht_keys;
static char *ht_used;

static u64 encode_sides(i64 a, i64 b, i64 c) {
    return ((u64)a << 42) | ((u64)b << 21) | (u64)c;
}

static unsigned hash_u64(u64 key) {
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return (unsigned)(key & HT_MASK);
}

static int seen_insert(u64 key) {
    unsigned h = hash_u64(key);
    while (ht_used[h]) {
        if (ht_keys[h] == key) return 0;  // already present
        h = (h + 1) & HT_MASK;
    }
    ht_used[h] = 1;
    ht_keys[h] = key;
    return 1;
}

long long p998_native(void) {
    build_partners(LIMIT);

    ht_keys = calloc(HT_SIZE, sizeof(u64));
    ht_used = calloc(HT_SIZE, 1);

    i64 total = 0;

    for (int m = 1; m <= LIMIT; m++) {
        i64 mm = (i64)m * m;
        int np = partner_count[m];
        i64 *plist = partner_ptr[m];

        // Build row: (0, m) followed by partners with their hypotenuses.
        // row[i] = (x, hx) where hx = isqrt(mm + x^2).
        // We need to sort by x. (0, m) is first, partners are already sorted.
        // But we need to merge (0, m) with the sorted partners.
        // Since partners are sorted and 0 <= all partners, (0,m) is first.

        int row_len = np + 1;
        i64 *row_x = malloc(row_len * sizeof(i64));
        i64 *row_h = malloc(row_len * sizeof(i64));
        row_x[0] = 0; row_h[0] = m;
        for (int i = 0; i < np; i++) {
            row_x[i + 1] = plist[i];
            row_h[i + 1] = isqrt_i64(mm + plist[i] * plist[i]);
        }

        // Edge-aligned minima.
        for (int i = 0; i < row_len; i++) {
            for (int j = i; j < row_len; j++) {
                i64 x = row_x[i], hx = row_h[i];
                i64 y = row_x[j], hy = row_h[j];
                i64 base = x + y;
                if (base == 0) continue;
                if (base > m) break;  // sorted, so no further j works
                if (x * y < m * (m - base)) continue;

                // Sort sides.
                i64 s[3] = {base, hx, hy};
                // Simple sort for 3 elements.
                if (s[0] > s[1]) { i64 t = s[0]; s[0] = s[1]; s[1] = t; }
                if (s[1] > s[2]) { i64 t = s[1]; s[1] = s[2]; s[2] = t; }
                if (s[0] > s[1]) { i64 t = s[0]; s[0] = s[1]; s[1] = t; }

                u64 key = encode_sides(s[0], s[1], s[2]);
                if (seen_insert(key)) {
                    total += s[0] + s[1] + s[2];
                }
            }
        }

        // Four-sided, balanced minima.
        for (int i = 0; i < row_len; i++) {
            for (int j = i; j < row_len; j++) {
                i64 u = row_x[i], hu = row_h[i];
                i64 v = row_x[j], hv = row_h[j];
                i64 twice_area = mm - u * v;
                if (twice_area <= 0) continue;
                i64 p = m - u;
                i64 q = m - v;
                i64 third2 = p * p + q * q;
                i64 third = isqrt_i64(third2);
                if (third * third != third2 || third == 0) continue;

                i64 s[3] = {third, hu, hv};
                if (s[0] > s[1]) { i64 t = s[0]; s[0] = s[1]; s[1] = t; }
                if (s[1] > s[2]) { i64 t = s[1]; s[1] = s[2]; s[2] = t; }
                if (s[0] > s[1]) { i64 t = s[0]; s[0] = s[1]; s[1] = t; }

                u64 key = encode_sides(s[0], s[1], s[2]);
                if (!seen_insert(key)) continue;

                if (is_minimum_square(s, twice_area, m)) {
                    total += s[0] + s[1] + s[2];
                }
            }
        }

        free(row_x);
        free(row_h);
    }

    // Cleanup.
    for (int m = 0; m <= LIMIT; m++) free(partner_ptr[m]);
    free(partner_ptr); free(partner_count); free(partner_cap);
    free(ht_keys); free(ht_used);

    return total;
}
