/* Project Euler 815: Group by Value
   E(60) = expected maximum number of non-empty piles.
   Port of the Python reference solver. */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Binomial coefficient C(n, k) for small n */
static long long comb(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    long long r = 1;
    for (int i = 0; i < k; i++) {
        r = r * (n - i) / (i + 1);
    }
    return r;
}

typedef struct {
    int n;
    long long total;
    int max_rem;
    double *inv;
    int *next0, *next1, *next2, *next3;
    unsigned short *num0, *num1, *num2, *num3;
    int start;
    /* sorted state indices by remaining */
    int *sorted_states;
    int *state_remaining;
    int *state_active;
} Model;

static void build_model(Model *m, int n) {
    m->n = n;
    m->total = comb(n + 4, 4);
    m->max_rem = 4 * n;

    m->inv = (double*)malloc(m->total * sizeof(double));
    m->next0 = (int*)calloc(m->total, sizeof(int));
    m->next1 = (int*)calloc(m->total, sizeof(int));
    m->next2 = (int*)calloc(m->total, sizeof(int));
    m->next3 = (int*)calloc(m->total, sizeof(int));
    m->num0 = (unsigned short*)calloc(m->total, sizeof(unsigned short));
    m->num1 = (unsigned short*)calloc(m->total, sizeof(unsigned short));
    m->num2 = (unsigned short*)calloc(m->total, sizeof(unsigned short));
    m->num3 = (unsigned short*)calloc(m->total, sizeof(unsigned short));

    /* prefix_s[s] = sum_{i=0..s-1} C(i+3, 3) */
    long long *prefix_s = (long long*)calloc(n + 2, sizeof(long long));
    for (int s = 0; s <= n; s++) {
        prefix_s[s + 1] = prefix_s[s] + comb(s + 3, 3);
    }

    /* x0_pref[s][x0] = prefix within fixed s for x0 */
    long long **x0_pref = (long long**)malloc((n + 1) * sizeof(long long*));
    for (int s = 0; s <= n; s++) {
        x0_pref[s] = (long long*)calloc(s + 2, sizeof(long long));
        long long running = 0;
        for (int x0 = 0; x0 <= s; x0++) {
            x0_pref[s][x0] = running;
            int t = s - x0;
            running += comb(t + 2, 2);
        }
        x0_pref[s][s + 1] = running;
    }

    /* idx function */
    /* base = prefix_s[s] + x0_pref[s][x0] + x1*(t+1) - x1*(x1-1)/2 + x2 */

    /* Enumerate all states in order, building transition data */
    int i = 0;
    m->state_remaining = (int*)malloc(m->total * sizeof(int));
    m->state_active = (int*)malloc(m->total * sizeof(int));

    for (int s = 0; s <= n; s++) {
        for (int x0 = 0; x0 <= s; x0++) {
            int t = s - x0;
            for (int x1 = 0; x1 <= t; x1++) {
                int u = t - x1;
                for (int x2 = 0; x2 <= u; x2++) {
                    int x3 = u - x2;
                    int remaining = 4 * x0 + 3 * x1 + 2 * x2 + x3;
                    int active = x1 + x2 + x3;

                    m->inv[i] = (remaining == 0) ? 0.0 : 1.0 / remaining;
                    m->state_remaining[i] = remaining;
                    m->state_active[i] = active;

                    /* Compute next state indices using idx function */
                    if (x0) {
                        /* next: (x0-1, x1+1, x2, x3) */
                        int ns = (x0 - 1) + (x1 + 1) + x2 + x3;
                        int nt = ns - (x0 - 1);
                        long long base = prefix_s[ns] + x0_pref[ns][x0 - 1];
                        base += (long long)(x1 + 1) * (nt + 1) - (long long)(x1 + 1) * x1 / 2;
                        base += x2;
                        m->next0[i] = (int)base;
                        m->num0[i] = 4 * x0;
                    }
                    if (x1) {
                        int ns = x0 + (x1 - 1) + (x2 + 1) + x3;
                        int nt = ns - x0;
                        long long base = prefix_s[ns] + x0_pref[ns][x0];
                        base += (long long)(x1 - 1) * (nt + 1) - (long long)(x1 - 1) * (x1 - 2) / 2;
                        base += (x2 + 1);
                        m->next1[i] = (int)base;
                        m->num1[i] = 3 * x1;
                    }
                    if (x2) {
                        int ns = x0 + x1 + (x2 - 1) + (x3 + 1);
                        int nt = ns - x0;
                        long long base = prefix_s[ns] + x0_pref[ns][x0];
                        base += (long long)x1 * (nt + 1) - (long long)x1 * (x1 - 1) / 2;
                        base += (x2 - 1);
                        m->next2[i] = (int)base;
                        m->num2[i] = 2 * x2;
                    }
                    if (x3) {
                        int ns = x0 + x1 + x2 + (x3 - 1);
                        int nt = ns - x0;
                        long long base = prefix_s[ns] + x0_pref[ns][x0];
                        base += (long long)x1 * (nt + 1) - (long long)x1 * (x1 - 1) / 2;
                        base += x2;
                        m->next3[i] = (int)base;
                        m->num3[i] = x3;
                    }

                    i++;
                }
            }
        }
    }

    /* Start state: (n, 0, 0, 0) */
    {
        int s = n, x0 = n, x1 = 0, x2 = 0, x3 = 0;
        int t = s - x0;
        long long base = prefix_s[s] + x0_pref[s][x0];
        base += (long long)x1 * (t + 1) - (long long)x1 * (x1 - 1) / 2;
        base += x2;
        m->start = (int)base;
    }

    /* Build sorted state indices by remaining */
    m->sorted_states = (int*)malloc(m->total * sizeof(int));
    for (i = 0; i < m->total; i++) m->sorted_states[i] = i;

    /* Counting sort by remaining (0..max_rem) */
    int *count = (int*)calloc(m->max_rem + 2, sizeof(int));
    for (i = 0; i < m->total; i++) count[m->state_remaining[i] + 1]++;
    for (int r = 0; r <= m->max_rem; r++) count[r + 1] += count[r];
    int *sorted = (int*)malloc(m->total * sizeof(int));
    for (i = 0; i < m->total; i++) {
        sorted[count[m->state_remaining[i]]++] = i;
    }
    free(m->sorted_states);
    m->sorted_states = sorted;
    free(count);

    /* Cleanup */
    for (int s = 0; s <= n; s++) free(x0_pref[s]);
    free(x0_pref);
    free(prefix_s);
}

static double probability_max_less_than_k(Model *m, int k, double *dp) {
    int total = m->total;
    int max_rem = m->max_rem;
    int max_active = k - 1;
    if (max_active > m->n) max_active = m->n;

    /* Reset dp */
    memset(dp, 0, total * sizeof(double));
    dp[0] = 1.0; /* terminal state */

    /* Process states in order of increasing remaining */
    for (int idx = 0; idx < total; idx++) {
        int s = m->sorted_states[idx];
        int rem = m->state_remaining[s];
        if (rem == 0) continue; /* terminal state already set */
        int active = m->state_active[s];
        if (active > max_active) continue; /* probability 0 */

        double val = 0.0;
        if (m->num0[s]) val += m->num0[s] * dp[m->next0[s]];
        if (m->num1[s]) val += m->num1[s] * dp[m->next1[s]];
        if (m->num2[s]) val += m->num2[s] * dp[m->next2[s]];
        if (m->num3[s]) val += m->num3[s] * dp[m->next3[s]];
        dp[s] = val * m->inv[s];
    }

    return dp[m->start];
}

static double expected_max_piles(int n) {
    Model m;
    build_model(&m, n);

    double *dp = (double*)malloc(m.total * sizeof(double));
    double total = 0.0;

    for (int k = 1; k <= n; k++) {
        double p = probability_max_less_than_k(&m, k, dp);
        total += (1.0 - p);
    }

    free(dp);
    free(m.inv);
    free(m.next0); free(m.next1); free(m.next2); free(m.next3);
    free(m.num0); free(m.num1); free(m.num2); free(m.num3);
    free(m.sorted_states);
    free(m.state_remaining);
    free(m.state_active);

    return total;
}

double p815_native(void) {
    return expected_max_piles(60);
}
