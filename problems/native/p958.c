// Project Euler 958: Euclid's Labour
// Meet-in-the-middle search on the Stern-Brocot tree.
#include <stdint.h>
#include <stdlib.h>

typedef long long i64;
typedef __int128 i128;

// Global state shared between consider() and check().
static i64 g_n;
static int g_steps;
static int g_split_depth;
static i64 g_best_steps;
static i64 g_best_value;

// State for check(), set by consider() at the midpoint.
static i64 g_local_a, g_local_b, g_norm;

// Modular inverse via extended GCD. Returns 1 if inverse exists, 0 otherwise.
static int mod_inverse(i64 value, i64 modulus, i64 *result) {
    value = ((value % modulus) + modulus) % modulus;
    i128 old_r = modulus, r = value;
    i128 old_s = 0, s = 1;
    while (r != 0) {
        i128 q = old_r / r;
        i128 tmp;
        tmp = old_r - q * r; old_r = r; r = tmp;
        tmp = old_s - q * s; old_s = s; s = tmp;
    }
    if (old_r != 1) return 0;
    *result = (i64)(((old_s % modulus) + modulus) % modulus);
    return 1;
}

// Check a candidate at the midpoint.
static void check(i64 candidate_a, i64 candidate_b, int current_depth) {
    if (candidate_a < 0 || candidate_b < 0)
        return;
    if (candidate_a * candidate_a + candidate_b * candidate_b > g_norm)
        return;

    i64 x = candidate_a, y = candidate_b;
    i64 vx = g_local_a, vy = g_local_b;

    if (g_steps & 1) {
        if ((vx & 1) || (y & 1))
            return;
        x -= y / 2;
        vy += vx / 2;
        if ((vy & 1) || (x & 1))
            return;
        x /= 2;
        y /= 2;
        vx /= 2;
        vy /= 2;
    }

    if (x * vx + y * vy != g_n)
        return;

    int remaining_steps = g_steps - g_split_depth;
    int used_steps = 0;
    while (used_steps <= remaining_steps && x != 0 && y != 0) {
        if (x > y) {
            i64 t = x; x = y; y = t;
            t = vx; vx = vy; vy = t;
        }
        y -= x;
        vx += vy;
        used_steps++;
    }

    if (used_steps > remaining_steps)
        return;

    i64 residue = ((vx + vy - g_n) % g_n + g_n) % g_n;
    i64 inv_residue;
    if (!mod_inverse(residue, g_n, &inv_residue))
        return;

    i64 total_steps = current_depth + used_steps;
    i64 n = g_n;
    i64 value = residue;
    if (n - residue < value) value = n - residue;
    if (inv_residue < value) value = inv_residue;
    if (n - inv_residue < value) value = n - inv_residue;

    if (total_steps < g_best_steps ||
        (total_steps == g_best_steps && value < g_best_value)) {
        g_best_steps = total_steps;
        g_best_value = value;
    }
}

// Recursive search on the Stern-Brocot tree.
static void consider(i64 basis_a, i64 basis_b, i64 coeff_a, i64 coeff_b,
                     int current_depth) {
    // Normalize: basis_a <= basis_b
    if (basis_a > basis_b) {
        i64 t = basis_a; basis_a = basis_b; basis_b = t;
        t = coeff_a; coeff_a = coeff_b; coeff_b = t;
    }

    // Adjust coefficients to be non-negative
    if (coeff_b < 0)
        return;
    if (coeff_a < 0) {
        i64 shift = (-coeff_a + basis_b - 1) / basis_b;
        coeff_a += shift * basis_b;
        coeff_b -= shift * basis_a;
        if (coeff_b < 0)
            return;
    }

    // Check dot product
    if (basis_a * coeff_a + basis_b * coeff_b != g_n)
        return;

    // Midpoint check
    if (current_depth == g_split_depth) {
        i64 local_a = basis_a, local_b = basis_b;
        i64 local_ca = coeff_a, local_cb = coeff_b;

        if (g_steps & 1) {
            local_a *= 2; local_b *= 2;
            local_ca *= 2; local_cb *= 2;
            local_b -= local_a / 2;
            local_ca += local_cb / 2;
            if (local_a * local_ca + local_b * local_cb != 4 * g_n)
                return;
        }

        i64 norm = local_a * local_a + local_b * local_b;
        if (norm < g_n)
            return;

        i64 cross = local_ca * local_b - local_cb * local_a;
        i64 shift = cross / norm;
        local_ca -= shift * local_b;
        local_cb += shift * local_a;
        cross -= shift * norm;

        if (cross < 0) {
            local_ca += local_b;
            local_cb -= local_a;
        }

        g_local_a = local_a;
        g_local_b = local_b;
        g_norm = norm;

        check(local_ca, local_cb, current_depth);
        check(local_ca - local_b, local_cb + local_a, current_depth);
        return;
    }

    // Pruning: compute basis at bottom of tree
    i64 x = basis_a, y = basis_b;
    for (int i = current_depth; i < g_steps / 2; i++) {
        if (x > y) { i64 t = x; x = y; y = t; }
        x += y;
        i64 t = x; x = y; y = t;
    }
    if (x > y) { i64 t = x; x = y; y = t; }

    if (g_steps & 1) {
        if (5 * y * y / 4 + x * y + x * x < g_n)
            return;
    } else {
        if (x * x + y * y < g_n)
            return;
    }

    // Recursive calls
    consider(basis_b, basis_a + basis_b, coeff_b - coeff_a, coeff_a,
             current_depth + 1);
    if (basis_a > 0 && basis_a < basis_b) {
        consider(basis_a, basis_a + basis_b, coeff_a - coeff_b, coeff_b,
                 current_depth + 1);
    }
}

static i64 f(i64 n) {
    g_n = n;
    int steps = 0;
    while (1) {
        g_steps = steps;
        g_split_depth = (steps + 1) / 2;
        g_best_steps = steps;
        g_best_value = n + 1;

        consider(0, 1, 0, n, 0);

        if (g_best_value <= n)
            return g_best_value;

        steps++;
    }
}

long long p958_native(void) {
    return f(1000000000000LL + 39);
}
