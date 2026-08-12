// Project Euler 882: Removing Bits
// Compute S(100000) for a combinatorial game where each number x has a game value
// g(x) that is a dyadic rational. We compute exactly using scaled integers.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// _simplest_between: return the simplest dyadic strictly between bounds.
// Values are integers scaled by DEN = 2^max_exp.
// a_scaled = lower bound (max of Left options), b_scaled = upper bound (min of Right options).
// -1 represents None (infinity).
static long long simplest_between(long long a_scaled, long long b_scaled, int max_exp) {
    long long den = 1LL << max_exp;

    if (a_scaled == -1) { // (-inf, b)
        if (b_scaled == -1)
            return 0;
        return ((b_scaled - 1) / den) * den;
    }

    if (b_scaled == -1) { // (a, +inf)
        return ((a_scaled / den) + 1) * den;
    }

    // a < b guaranteed by caller
    for (int k = 0; k <= max_exp; k++) {
        long long step = 1LL << (max_exp - k);
        long long m_low = a_scaled / step + 1;
        long long m_high = (b_scaled - 1) / step;
        if (m_low <= m_high)
            return m_low * step;
    }
    // Should never reach here
    return 0;
}

static long long compute_S(long long n) {
    if (n <= 1) return n;

    int max_exp = 0;
    long long tmp = n;
    while (tmp > 0) { max_exp++; tmp >>= 1; }
    long long den = 1LL << max_exp;

    long long *g = calloc(n + 1, sizeof(long long));
    long long total_scaled = 0;

    // Temporary arrays for options
    long long *left_opts = malloc(64 * sizeof(long long));
    long long *right_opts = malloc(64 * sizeof(long long));

    for (long long x = 1; x <= n; x++) {
        int left_count = 0, right_count = 0;

        int bits = 0;
        long long tmpx = x;
        while (tmpx > 0) { bits++; tmpx >>= 1; }

        for (int j = 0; j < bits; j++) {
            int bit = (x >> j) & 1;
            long long higher = x >> (j + 1);
            long long lower = x & ((1LL << j) - 1);
            long long y = (higher << j) | lower;

            if (bit) {
                left_opts[left_count++] = g[y];
            } else {
                right_opts[right_count++] = g[y];
            }
        }

        // Canonical pruning
        int changed = 1;
        while (changed) {
            changed = 0;
            if (right_count > 0) {
                long long min_r = right_opts[0];
                for (int i = 1; i < right_count; i++)
                    if (right_opts[i] < min_r) min_r = right_opts[i];
                int new_left_count = 0;
                for (int i = 0; i < left_count; i++)
                    if (left_opts[i] < min_r)
                        left_opts[new_left_count++] = left_opts[i];
                if (new_left_count != left_count) {
                    left_count = new_left_count;
                    changed = 1;
                }
            }
            if (left_count > 0) {
                long long max_l = left_opts[0];
                for (int i = 1; i < left_count; i++)
                    if (left_opts[i] > max_l) max_l = left_opts[i];
                int new_right_count = 0;
                for (int i = 0; i < right_count; i++)
                    if (right_opts[i] > max_l)
                        right_opts[new_right_count++] = right_opts[i];
                if (new_right_count != right_count) {
                    right_count = new_right_count;
                    changed = 1;
                }
            }
        }

        long long max_l = -1, min_r = -1;
        if (left_count > 0) {
            max_l = left_opts[0];
            for (int i = 1; i < left_count; i++)
                if (left_opts[i] > max_l) max_l = left_opts[i];
        }
        if (right_count > 0) {
            min_r = right_opts[0];
            for (int i = 1; i < right_count; i++)
                if (right_opts[i] < min_r) min_r = right_opts[i];
        }

        g[x] = simplest_between(max_l, min_r, max_exp);
        total_scaled += x * g[x];
    }

    free(left_opts);
    free(right_opts);
    free(g);

    return (total_scaled + den - 1) / den;
}

long long p882_native(void) {
    return compute_S(100000);
}
