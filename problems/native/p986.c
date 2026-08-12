// Project Euler 986: Sum of G(c,d) for 1<=c,d<=160.
// Cellular automaton threshold simulation with cubic extrapolation.
// Port of the Python reference solver to C.
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef long long i64;

static int extinct_for_k1(int n, i64 k) {
    if (k == 0) return 1;

    int size = n + 1;
    int last = size - 1;
    i64 *cells = (i64 *)calloc((size_t)size, sizeof(i64));
    cells[last] = k;
    int zero_count = last;

    while (1) {
        for (int i = 0; i < last; i++) {
            i64 old = cells[i];
            i64 nxt = (old + cells[i + 1]) >> 1;
            cells[i] = nxt;
            if (old) {
                if (!nxt) zero_count++;
            } else if (nxt) {
                zero_count--;
            }
        }

        i64 old = cells[last];
        i64 nxt = (old + cells[0]) >> 1;
        cells[last] = nxt;
        if (old) {
            if (!nxt) zero_count++;
        } else if (nxt) {
            zero_count--;
        }

        if (zero_count == size) {
            free(cells);
            return 1;
        }
        if (zero_count == 0) {
            free(cells);
            return 0;
        }
    }
}

static i64 threshold_k1_plain(int n) {
    i64 lo = 0, hi = 1;
    while (extinct_for_k1(n, hi)) {
        lo = hi;
        hi *= 2;
    }
    while (lo + 1 < hi) {
        i64 mid = (lo + hi) / 2;
        if (extinct_for_k1(n, mid))
            lo = mid;
        else
            hi = mid;
    }
    return lo;
}

static i64 predict_k1_from_previous(i64 *s, int n) {
    i64 a = s[n - 32];
    i64 b = s[n - 24];
    i64 c = s[n - 16];
    i64 d = s[n - 8];
    return d + (d - c) + (d - 2 * c + b) + (d - 3 * c + 3 * b - a);
}

static i64 threshold_k1_with_guess(int n, i64 guess) {
    i64 lo = guess - 4096;
    if (lo < 0) lo = 0;
    i64 hi = guess + 4096;

    while (lo > 0 && !extinct_for_k1(n, lo)) {
        hi = lo;
        lo /= 2;
    }
    while (extinct_for_k1(n, hi)) {
        lo = hi;
        hi *= 2;
    }
    while (lo + 1 < hi) {
        i64 mid = (lo + hi) / 2;
        if (extinct_for_k1(n, mid))
            lo = mid;
        else
            hi = mid;
    }
    return lo;
}

static i64 gcd_ll(i64 a, i64 b) {
    while (b) {
        i64 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

/* Exception table for d=1: H(c,1) exceptions */
static const int exc_c[] = {2, 3, 4, 5, 6, 8, 10};
static const int exc_h[] = {3, 5, 7, 11, 13, 21, 31};
static const int exc_count = 7;

static int get_exception(int c) {
    for (int i = 0; i < exc_count; i++) {
        if (exc_c[i] == c) return exc_h[i];
    }
    return -1;
}

long long p986_native(void) {
    int limit = 160;
    int max_n = limit + (limit - 1) / 2; /* 239 */

    i64 *s = (i64 *)calloc((size_t)(max_n + 1), sizeof(i64));

    for (int n = 1; n <= max_n; n++) {
        if (n < 33) {
            s[n] = threshold_k1_plain(n);
        } else {
            i64 guess = predict_k1_from_previous(s, n);
            s[n] = threshold_k1_with_guess(n, guess);
        }
    }

    /* Memoization table for reduced pairs (cr, dr), both <= limit */
    int dim = limit + 1;
    i64 *memo = (i64 *)malloc((size_t)dim * dim * sizeof(i64));
    memset(memo, -1, (size_t)dim * dim * sizeof(i64));

    i64 total = 0;
    for (int c = 1; c <= limit; c++) {
        for (int d = 1; d <= limit; d++) {
            i64 g = gcd_ll(c, d);
            int cr = (int)(c / g);
            int dr = (int)(d / g);
            i64 *slot = &memo[cr * dim + dr];
            i64 val = *slot;
            if (val < 0) {
                i64 h;
                if (dr == 1) {
                    int exc = get_exception(cr);
                    if (exc >= 0)
                        h = exc;
                    else
                        h = s[dr + (cr - 1) / 2];
                } else {
                    h = s[dr + (cr - 1) / 2];
                }
                val = 2 * h + 1;
                *slot = val;
            }
            total += val;
        }
    }

    free(s);
    free(memo);
    return total;
}
