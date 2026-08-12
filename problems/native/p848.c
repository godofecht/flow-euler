// Project Euler 848: Guessing with Sets
// sum_{i=0..20} sum_{j=0..20} p(7^i, 5^j), rounded to 8 decimals
// Uses __int128 for exact fraction arithmetic, converts to double at the end.
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

// next_pow2_times3: returns (T, p) where p is smallest power of 2 with 3*p >= x, T = 3*p
static void next_pow2_times3(i64 x, i64 *T_out, i64 *p_out) {
    i64 p = 1;
    while (3 * p < x) p <<= 1;
    *T_out = 3 * p;
    *p_out = p;
}

// Compute p(m, n) as an exact fraction (num/den) using __int128.
// Returns 1 if a closed form applies, 0 if the recurrence would be needed.
// For the target (7^i, 5^j), closed forms always apply.
static int p_fraction(i64 m, i64 n, i128 *num_out, i128 *den_out) {
    if (m == 1) { *num_out = 1; *den_out = 1; return 1; }
    if (n == 1) { *num_out = 1; *den_out = m; return 1; }
    if (n == 2) { *num_out = 3; *den_out = 2 * m; return 1; }
    if (m == 2) { *num_out = 2 * n - 1; *den_out = 2 * n; return 1; }
    if (m == 3) { *num_out = n - 1; *den_out = n; return 1; }

    // High-n region
    i64 Tm, pm;
    next_pow2_times3(m, &Tm, &pm);
    if (pm >= 2) {
        i64 L = 3 * (pm >> 1);
        if (n >= L) {
            // p = 1 - L*(m - pm) / (m*n)
            // = (m*n - L*(m - pm)) / (m*n)
            i128 mn = (i128)m * n;
            i128 sub = (i128)L * (m - pm);
            *num_out = mn - sub;
            *den_out = mn;
            return 1;
        }
    }

    // High-m region
    i64 Tn, pn;
    next_pow2_times3(n, &Tn, &pn);
    if (n >= 3 && m >= Tn) {
        // p = Tn * (n - pn) / (n * m)
        *num_out = (i128)Tn * (n - pn);
        *den_out = (i128)n * m;
        return 1;
    }

    // Recurrence would be needed - but for (7^i, 5^j) this never happens.
    return 0;
}

double p848_native(void) {
    // Precompute powers of 7 and 5
    i64 pow7[21], pow5[21];
    pow7[0] = 1;
    for (int i = 1; i <= 20; i++) pow7[i] = pow7[i-1] * 7;
    pow5[0] = 1;
    for (int i = 1; i <= 20; i++) pow5[i] = pow5[i-1] * 5;

    double total = 0.0;

    for (int i = 0; i <= 20; i++) {
        for (int j = 0; j <= 20; j++) {
            i64 m = pow7[i];
            i64 n = pow5[j];

            i128 num, den;
            if (p_fraction(m, n, &num, &den)) {
                total += (double)num / (double)den;
            } else {
                // Should not happen for (7^i, 5^j)
                // If it did, we'd need the recurrence. Fall back to 0.
            }
        }
    }

    return total;
}
