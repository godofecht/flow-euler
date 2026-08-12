// Project Euler 852: Coins in a Box
// S(50) rounded to 6 decimals
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_FLIPS 180
#define N_MAX 50

static double stop_value(double p) {
    if (p >= 0.5) return 70.0 * p - 50.0;
    return 20.0 - 70.0 * p;
}

static double compute_g_for_fraction(int a, int b, const double *pow3, const double *inv2) {
    if (a <= 0 || a >= b) return 20.0;

    double p0 = (double)a / (double)b;
    double odds0 = p0 / (1.0 - p0);

    double *next_row = (double*)calloc(MAX_FLIPS + 2, sizeof(double));
    double *curr_row = (double*)calloc(MAX_FLIPS + 2, sizeof(double));

    // Initialize at n = MAX_FLIPS with forced stop
    double base = odds0 * inv2[MAX_FLIPS];
    for (int h = 0; h <= MAX_FLIPS; h++) {
        double odds = base * pow3[h];
        double p = odds / (1.0 + odds);
        next_row[h] = stop_value(p);
    }

    // Backward induction
    for (int n = MAX_FLIPS - 1; n >= 0; n--) {
        base = odds0 * inv2[n];
        for (int h = 0; h <= n; h++) {
            double odds = base * pow3[h];
            double p = odds / (1.0 + odds);
            double stop_v = stop_value(p);
            double qh = 0.5 + 0.25 * p;
            double cont_v = -1.0 + qh * next_row[h + 1] + (1.0 - qh) * next_row[h];
            curr_row[h] = (stop_v >= cont_v) ? stop_v : cont_v;
        }
        double *tmp = next_row;
        next_row = curr_row;
        curr_row = tmp;
    }

    double result = next_row[0];
    free(next_row);
    free(curr_row);
    return result;
}

// GCD
static int gcd_int(int a, int b) {
    while (b) { int t = a % b; a = b; b = t; }
    return a;
}

double p852_native(void) {
    int n = N_MAX;

    // Precompute pow3 and inv2
    double pow3[MAX_FLIPS + 1];
    double inv2[MAX_FLIPS + 1];
    pow3[0] = 1.0;
    inv2[0] = 1.0;
    for (int i = 1; i <= MAX_FLIPS; i++) {
        pow3[i] = pow3[i - 1] * 3.0;
        inv2[i] = inv2[i - 1] * 0.5;
    }

    // Collect all needed (a, b) fractions
    // For u in 0..n, f in 0..n, total = u+f, if total > 0: (u/g, total/g)
    // We need a cache. Use a hash map or just a 2D array keyed by (a, b) with b <= 2*n.
    // Max b = 2*n = 100. So we can use a 2D array.
    // a ranges 0..b, b ranges 1..2*n.
    int max_b = 2 * n;
    double **g_cache = (double**)malloc((size_t)(max_b + 1) * sizeof(double*));
    for (int b = 0; b <= max_b; b++)
        g_cache[b] = (double*)calloc((size_t)(b + 1), sizeof(double));

    // Track which (a,b) have been computed
    char **g_done = (char**)malloc((size_t)(max_b + 1) * sizeof(char*));
    for (int b = 0; b <= max_b; b++)
        g_done[b] = (char*)calloc((size_t)(b + 1), sizeof(char));

    for (int u = 0; u <= n; u++) {
        for (int f = 0; f <= n; f++) {
            int total = u + f;
            if (total == 0) continue;
            int g = gcd_int(u, total);
            int a = u / g;
            int b = total / g;
            if (!g_done[b][a]) {
                g_cache[b][a] = compute_g_for_fraction(a, b, pow3, inv2);
                g_done[b][a] = 1;
            }
        }
    }

    // V[u][f] = expected optimal score
    double **V = (double**)malloc((size_t)(n + 1) * sizeof(double*));
    for (int u = 0; u <= n; u++)
        V[u] = (double*)calloc((size_t)(n + 1), sizeof(double));

    // Fill by increasing total = u+f
    for (int total = 1; total <= 2 * n; total++) {
        int u_min = (total - n > 0) ? total - n : 0;
        int u_max = (n < total) ? n : total;
        for (int u = u_min; u <= u_max; u++) {
            int f = total - u;
            int g = gcd_int(u, total);
            double immediate = g_cache[total / g][u / g];
            double exp_next = 0.0;
            if (u) exp_next += ((double)u / (double)total) * V[u - 1][f];
            if (f) exp_next += ((double)f / (double)total) * V[u][f - 1];
            V[u][f] = immediate + exp_next;
        }
    }

    double ans = V[n][n];

    // Cleanup
    for (int b = 0; b <= max_b; b++) { free(g_cache[b]); free(g_done[b]); }
    free(g_cache);
    free(g_done);
    for (int u = 0; u <= n; u++) free(V[u]);
    free(V);

    return ans;
}
