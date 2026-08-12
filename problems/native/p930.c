/* Project Euler 930: The Gathering
 * Compute G(12, 12) = sum_{n=2..12} sum_{m=2..12} F(n, m)
 * where F(n, m) = sum_{k != 0} 1/(1 - lambda_k)
 * and lambda_k are eigenvalues of the random walk on (Z_n)^(m-1).
 *
 * The eigenvalues are parametrized by frequency vectors k in (Z_n)^(m-1).
 * We enumerate compositions of d = m-1 into n parts (counts of each residue),
 * computing sum_cos and sum_mod, with multinomial multiplicity.
 *
 * Output in scientific notation with 12 decimal places, no '+' in exponent.
 */

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

static int BINOM[12][12];

static void init_binom(void) {
    for (int n = 0; n < 12; n++) {
        BINOM[n][0] = 1;
        for (int k = 1; k <= n; k++)
            BINOM[n][k] = BINOM[n-1][k-1] + (k <= n-1 ? BINOM[n-1][k] : 0);
    }
}

/* Kahan summation state */
typedef struct { double total; double c; } Kahan;

static void kahan_add(Kahan *k, double value) {
    double y = value - k->c;
    double t = k->total + y;
    k->c = (t - k->total) - y;
    k->total = t;
}

/* Recursion state for F(n, m) */
typedef struct {
    int n, d;
    double cos_table[16];
    Kahan kah;
} FState;

static void rec(FState *st, int r, int remaining, int mult,
                double sum_cos, int sum_mod, int any_nonzero) {
    if (r == st->n - 1) {
        int cnt = remaining;
        double sum_cos2 = sum_cos + cnt * st->cos_table[r];
        int sum_mod2 = (sum_mod + cnt * r) % st->n;
        int any2 = any_nonzero || (cnt > 0 && r != 0);
        if (!any2) return; /* all-zero frequency vector, eigenvalue = 1 */
        double lam = (sum_cos2 + st->cos_table[sum_mod2]) / (double)(st->d + 1);
        double term = (double)mult / (1.0 - lam);
        kahan_add(&st->kah, term);
        return;
    }

    int *row = BINOM[remaining];
    double cr = st->cos_table[r];
    for (int cnt = 0; cnt <= remaining; cnt++) {
        rec(st, r + 1, remaining - cnt, mult * row[cnt],
            sum_cos + cnt * cr,
            (sum_mod + cnt * r) % st->n,
            any_nonzero || (cnt > 0 && r != 0));
    }
}

static double F(int n, int m) {
    int d = m - 1;
    FState st;
    st.n = n;
    st.d = d;
    st.kah.total = 0.0;
    st.kah.c = 0.0;

    for (int r = 0; r < n; r++)
        st.cos_table[r] = cos(2.0 * M_PI * r / n);

    rec(&st, 0, d, 1, 0.0, 0, 0);
    return st.kah.total;
}

static double G(int N, int M) {
    Kahan kah = {0.0, 0.0};
    for (int n = 2; n <= N; n++)
        for (int m = 2; m <= M; m++)
            kahan_add(&kah, F(n, m));
    return kah.total;
}

int p930_native(void) {
    init_binom();
    double ans = G(12, 12);

    /* Format as %.12e and strip '+' from exponent */
    char buf[64];
    snprintf(buf, sizeof(buf), "%.12e", ans);
    /* Find 'e' and strip '+' if present */
    char *e = strchr(buf, 'e');
    if (e) {
        char *plus = strchr(e, '+');
        if (plus) {
            /* Shift remaining chars left by 1 to remove '+' */
            memmove(plus, plus + 1, strlen(plus + 1) + 1);
        }
    }
    printf("%s\n", buf);
    return 0;
}
