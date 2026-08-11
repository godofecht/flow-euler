/*
 * Project Euler 794 - Seventeen Points
 *
 * DFS over all feasible insertion histories up to n=18.
 * All boundaries k/n are represented as integers scaled by D = lcm(1..18).
 * The minimal sum at n=17 is recorded, then divided by D and rounded
 * to 12 decimal places (ROUND_HALF_UP).
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>

static long long gcd_ll(long long a, long long b) {
    while (b) { long long t = b; b = a % b; a = t; }
    return a;
}

static long long lcm_ll(long long a, long long b) {
    return a / gcd_ll(a, b) * b;
}

#define MAX_N 18
#define TARGET_N 17

static long long g_scale[MAX_N + 1];
static long long g_best;
static int g_best_found;
static int g_exists_max;

static void rec(int n, int *order, long long *L, long long *U, long long sumL) {
    if (n == TARGET_N) {
        if (!g_best_found || sumL < g_best) {
            g_best = sumL;
            g_best_found = 1;
        }
    }
    if (n == MAX_N) {
        g_exists_max = 1;
        return;
    }

    int m = n + 1;
    long long sc = g_scale[m];

    int new_order[MAX_N + 1];
    long long newL[MAX_N + 1];
    long long newU[MAX_N + 1];

    for (int pos = 0; pos < m; pos++) {
        for (int i = 0; i < pos; i++) new_order[i] = order[i];
        new_order[pos] = m;
        for (int i = pos; i < n; i++) new_order[i + 1] = order[i];

        memcpy(newL, L, sizeof(long long) * (MAX_N + 1));
        memcpy(newU, U, sizeof(long long) * (MAX_N + 1));
        long long new_sumL = sumL;
        int feasible = 1;

        for (int k = 0; k < m; k++) {
            int pid = new_order[k];
            long long lb = (long long)k * sc;
            long long ub = (long long)(k + 1) * sc;

            if (lb > newL[pid]) {
                new_sumL += lb - newL[pid];
                newL[pid] = lb;
            }
            if (ub < newU[pid]) {
                newU[pid] = ub;
            }
            if (newL[pid] >= newU[pid]) {
                feasible = 0;
                break;
            }
        }

        if (feasible) {
            rec(m, new_order, newL, newU, new_sumL);
        }
    }
}

double p794_native(void) {
    long long denom = 1;
    for (int i = 1; i <= MAX_N; i++) {
        denom = lcm_ll(denom, i);
    }
    for (int n = 1; n <= MAX_N; n++) {
        g_scale[n] = denom / n;
    }

    long long L[MAX_N + 1];
    long long U[MAX_N + 1];
    for (int i = 0; i <= MAX_N; i++) {
        L[i] = 0;
        U[i] = denom;
    }

    int order[MAX_N + 1];
    order[0] = 1;

    g_best_found = 0;
    g_best = 0;
    g_exists_max = 0;

    rec(1, order, L, U, 0);

    /* Round g_best / denom to 12 decimal places, ROUND_HALF_UP. */
    __int128 num = (__int128)g_best * 1000000000000LL;
    __int128 q = num / denom;
    __int128 r = num % denom;
    if (2 * r >= (__int128)denom) q++;

    return (double)(long long)q / 1000000000000.0;
}
