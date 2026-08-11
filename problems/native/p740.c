// Project Euler 740: Secret Santa
// Probability that the last person ends up with at least one of their own slips.
// Port of the Python reference solver to C.
//
// State: (u1, u2, k) where
//   u1 = unprocessed non-last people with 1 slip in the hat
//   u2 = unprocessed non-last people with 2 slips in the hat
//   k  = slips of the last person remaining in the hat (0..2)
// sp (processed slips) is derived: sp = T_total - u1 - 2*u2 - k.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 100
#define MAXU (N + 2)
#define MAXK 3
#define STATE_SIZE (MAXU * MAXU * MAXK)

static inline int sidx(int u1, int u2, int k) {
    return (u1 * MAXU + u2) * MAXK + k;
}

static double q(int n) {
    if (n < 2) return 0.0;
    if (n == 2) return 1.0;

    double *cur = calloc(STATE_SIZE, sizeof(double));
    double *nxt = calloc(STATE_SIZE, sizeof(double));
    if (!cur || !nxt) { fprintf(stderr, "oom\n"); exit(1); }

    cur[sidx(0, n - 1, 2)] = 1.0;

    for (int t = 0; t < n - 1; t++) {
        int m = (n - 1) - t;
        int T_total = 2 * n - 2 * t;

        memset(nxt, 0, STATE_SIZE * sizeof(double));

        for (int u1 = 0; u1 <= m; u1++) {
            for (int u2 = 0; u2 <= m - u1; u2++) {
                int u0 = m - u1 - u2;
                for (int k = 0; k <= 2; k++) {
                    double prob = cur[sidx(u1, u2, k)];
                    if (prob == 0.0) continue;

                    int sp = T_total - u1 - 2 * u2 - k;
                    if (sp < 0) continue;

                    for (int s = 0; s <= 2; s++) {
                        int cnt = (s == 0) ? u0 : (s == 1) ? u1 : u2;
                        if (cnt == 0) continue;
                        double p_actor = prob * ((double)cnt / (double)m);

                        int uu1 = u1, uu2 = u2;
                        if (s == 1) uu1 -= 1;
                        else if (s == 2) uu2 -= 1;

                        int C1 = T_total - s;
                        double invC1 = 1.0 / (double)C1;
                        int C2 = (T_total - 1) - s;
                        double invC2 = 1.0 / (double)C2;

                        /* First draw: last-person slip */
                        if (k > 0) {
                            double p1 = (double)k * invC1;
                            int a = uu1, b = uu2, kk = k - 1, ss = sp;
                            if (kk > 0)
                                nxt[sidx(a, b, kk - 1)] += p_actor * p1 * ((double)kk * invC2);
                            if (a > 0)
                                nxt[sidx(a - 1, b, kk)] += p_actor * p1 * ((double)a * invC2);
                            if (b > 0)
                                nxt[sidx(a + 1, b - 1, kk)] += p_actor * p1 * (2.0 * (double)b * invC2);
                            if (ss > 0)
                                nxt[sidx(a, b, kk)] += p_actor * p1 * ((double)ss * invC2);
                        }
                        /* First draw: unprocessed(1) */
                        if (uu1 > 0) {
                            double p1 = (double)uu1 * invC1;
                            int a = uu1 - 1, b = uu2, kk = k, ss = sp;
                            if (kk > 0)
                                nxt[sidx(a, b, kk - 1)] += p_actor * p1 * ((double)kk * invC2);
                            if (a > 0)
                                nxt[sidx(a - 1, b, kk)] += p_actor * p1 * ((double)a * invC2);
                            if (b > 0)
                                nxt[sidx(a + 1, b - 1, kk)] += p_actor * p1 * (2.0 * (double)b * invC2);
                            if (ss > 0)
                                nxt[sidx(a, b, kk)] += p_actor * p1 * ((double)ss * invC2);
                        }
                        /* First draw: unprocessed(2) */
                        if (uu2 > 0) {
                            double p1 = 2.0 * (double)uu2 * invC1;
                            int a = uu1 + 1, b = uu2 - 1, kk = k, ss = sp;
                            if (kk > 0)
                                nxt[sidx(a, b, kk - 1)] += p_actor * p1 * ((double)kk * invC2);
                            if (a > 0)
                                nxt[sidx(a - 1, b, kk)] += p_actor * p1 * ((double)a * invC2);
                            if (b > 0)
                                nxt[sidx(a + 1, b - 1, kk)] += p_actor * p1 * (2.0 * (double)b * invC2);
                            if (ss > 0)
                                nxt[sidx(a, b, kk)] += p_actor * p1 * ((double)ss * invC2);
                        }
                        /* First draw: processed pool */
                        if (sp > 0) {
                            double p1 = (double)sp * invC1;
                            int a = uu1, b = uu2, kk = k, ss = sp - 1;
                            if (kk > 0)
                                nxt[sidx(a, b, kk - 1)] += p_actor * p1 * ((double)kk * invC2);
                            if (a > 0)
                                nxt[sidx(a - 1, b, kk)] += p_actor * p1 * ((double)a * invC2);
                            if (b > 0)
                                nxt[sidx(a + 1, b - 1, kk)] += p_actor * p1 * (2.0 * (double)b * invC2);
                            if (ss > 0)
                                nxt[sidx(a, b, kk)] += p_actor * p1 * ((double)ss * invC2);
                        }
                    }
                }
            }
        }

        double *tmp = cur; cur = nxt; nxt = tmp;
    }

    /* After n-1 people, only 2 slips remain. Fail iff k > 0. */
    double ans = 0.0;
    for (int i = 0; i < STATE_SIZE; i++) {
        int k = i % MAXK;
        if (k > 0) ans += cur[i];
    }

    free(cur);
    free(nxt);
    return ans;
}

double p740_native(void) {
    return q(100);
}
