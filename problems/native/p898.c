// Project Euler 898: Claire Voyant
// Optimal probability that Claire guesses the coin toss.
// Meet-in-the-middle with exact rational likelihood ratios via GMP.
// 51 students with lie probabilities 25%..75%, paired into 25 3-outcome variables.
#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int n_outcomes;
    int mul_n[3], mul_d[3];
    double ph[3], pt[3];
} Variable;

static Variable variables[60];
static int n_vars;

static void build_variables(const int *ps, int n_ps) {
    int counts[101];
    memset(counts, 0, sizeof(counts));
    for (int i = 0; i < n_ps; i++) counts[ps[i]]++;

    n_vars = 0;
    for (int k = 0; k < 50; k++) {
        int m = counts[k] < counts[100 - k] ? counts[k] : counts[100 - k];
        if (m <= 0) continue;
        counts[k] -= m;
        counts[100 - k] -= m;

        int a = 100 - k, b = k;
        int num_plus = a * a, den_plus = b * b;
        int num_minus = b * b, den_minus = a * a;

        double p_plus_H = (double)(a * a) / 10000.0;
        double p_zero = (double)(2 * a * b) / 10000.0;
        double p_minus_H = (double)(b * b) / 10000.0;

        for (int r = 0; r < m; r++) {
            Variable *v = &variables[n_vars++];
            v->n_outcomes = 3;
            v->mul_n[0] = num_minus; v->mul_d[0] = den_minus;
            v->ph[0] = p_minus_H; v->pt[0] = p_plus_H;
            v->mul_n[1] = 1; v->mul_d[1] = 1;
            v->ph[1] = p_zero; v->pt[1] = p_zero;
            v->mul_n[2] = num_plus; v->mul_d[2] = den_plus;
            v->ph[2] = p_plus_H; v->pt[2] = p_minus_H;
        }
    }

    counts[50] = 0;

    for (int k = 0; k <= 100; k++) {
        for (int r = 0; r < counts[k]; r++) {
            int a = 100 - k, b = k;
            Variable *v = &variables[n_vars++];
            v->n_outcomes = 2;
            v->mul_n[0] = b; v->mul_d[0] = a;
            v->ph[0] = b / 100.0; v->pt[0] = a / 100.0;
            v->mul_n[1] = a; v->mul_d[1] = b;
            v->ph[1] = a / 100.0; v->pt[1] = b / 100.0;
        }
    }
}

#define KEY_BYTES 64

typedef struct {
    unsigned char key[KEY_BYTES];
    double ph, pt;
} RState;

static RState *r_states;
static int n_r;
static double *suf_h, *suf_t;

static void mpz_to_bytes(unsigned char *buf, const mpz_t val) {
    memset(buf, 0, KEY_BYTES);
    size_t count;
    mpz_export(buf, &count, 1, 1, 1, 0, val);
    if (count > 0 && (int)count < KEY_BYTES) {
        memmove(buf + (KEY_BYTES - count), buf, count);
        memset(buf, 0, KEY_BYTES - count);
    }
}

static int rstate_cmp(const void *a, const void *b) {
    return memcmp(((const RState *)a)->key, ((const RState *)b)->key, KEY_BYTES);
}

static int bisect_left(const unsigned char *threshold) {
    int lo = 0, hi = n_r;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (memcmp(r_states[mid].key, threshold, KEY_BYTES) < 0)
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo;
}

static void enumerate_right(Variable *vars, int n, int shift) {
    mpz_t *nums, *dens;
    double *phs, *pts;
    int count = 1;
    nums = malloc(sizeof(mpz_t)); mpz_init_set_ui(nums[0], 1);
    dens = malloc(sizeof(mpz_t)); mpz_init_set_ui(dens[0], 1);
    phs = malloc(sizeof(double)); phs[0] = 1.0;
    pts = malloc(sizeof(double)); pts[0] = 1.0;

    for (int v = 0; v < n; v++) {
        int no = vars[v].n_outcomes;
        int new_count = count * no;
        mpz_t *nn = malloc(new_count * sizeof(mpz_t));
        mpz_t *nd = malloc(new_count * sizeof(mpz_t));
        double *nph = malloc(new_count * sizeof(double));
        double *npt = malloc(new_count * sizeof(double));

        int idx = 0;
        for (int i = 0; i < count; i++) {
            for (int o = 0; o < no; o++) {
                mpz_init(nn[idx]);
                mpz_mul_ui(nn[idx], nums[i], vars[v].mul_n[o]);
                mpz_init(nd[idx]);
                mpz_mul_ui(nd[idx], dens[i], vars[v].mul_d[o]);
                nph[idx] = phs[i] * vars[v].ph[o];
                npt[idx] = pts[i] * vars[v].pt[o];
                idx++;
            }
        }

        for (int i = 0; i < count; i++) { mpz_clear(nums[i]); mpz_clear(dens[i]); }
        free(nums); free(dens); free(phs); free(pts);
        nums = nn; dens = nd; phs = nph; pts = npt;
        count = new_count;
    }

    n_r = count;
    r_states = malloc(count * sizeof(RState));
    mpz_t shifted;
    mpz_init(shifted);
    for (int i = 0; i < count; i++) {
        mpz_mul_2exp(shifted, nums[i], shift);
        mpz_fdiv_q(shifted, shifted, dens[i]);
        mpz_to_bytes(r_states[i].key, shifted);
        r_states[i].ph = phs[i];
        r_states[i].pt = pts[i];
        mpz_clear(nums[i]); mpz_clear(dens[i]);
    }
    mpz_clear(shifted);
    free(nums); free(dens); free(phs); free(pts);

    qsort(r_states, n_r, sizeof(RState), rstate_cmp);

    suf_h = malloc((n_r + 1) * sizeof(double));
    suf_t = malloc((n_r + 1) * sizeof(double));
    suf_h[n_r] = 0.0;
    suf_t[n_r] = 0.0;
    for (int i = n_r - 1; i >= 0; i--) {
        suf_h[i] = suf_h[i + 1] + r_states[i].ph;
        suf_t[i] = suf_t[i + 1] + r_states[i].pt;
    }
}

static double pHA, cH, pTA, cT;
static int n_left, shift_g;
static mpz_t lnum[14], lden[14];

static void enum_left(int vi, double ph, double pt) {
    if (vi == n_left) {
        mpz_t shifted, threshold;
        mpz_init(shifted); mpz_init(threshold);
        if (n_left == 0) {
            mpz_set_ui(threshold, 1);
        } else {
            mpz_mul_2exp(shifted, lden[vi - 1], shift_g);
            mpz_fdiv_q(threshold, shifted, lnum[vi - 1]);
        }
        unsigned char thr_bytes[KEY_BYTES];
        mpz_to_bytes(thr_bytes, threshold);
        int j = bisect_left(thr_bytes);

        double y, t;
        y = ph * suf_h[j] - cH;
        t = pHA + y; cH = (t - pHA) - y; pHA = t;
        y = pt * suf_t[j] - cT;
        t = pTA + y; cT = (t - pTA) - y; pTA = t;

        mpz_clear(shifted); mpz_clear(threshold);
        return;
    }

    Variable *v = &variables[vi];
    for (int o = 0; o < v->n_outcomes; o++) {
        if (vi == 0) {
            mpz_set_ui(lnum[0], v->mul_n[o]);
            mpz_set_ui(lden[0], v->mul_d[o]);
        } else {
            mpz_mul_ui(lnum[vi], lnum[vi - 1], v->mul_n[o]);
            mpz_mul_ui(lden[vi], lden[vi - 1], v->mul_d[o]);
        }
        enum_left(vi + 1, ph * v->ph[o], pt * v->pt[o]);
    }
}

static void compute_half_bound(mpz_t result, int start, int end) {
    mpz_set_ui(result, 1);
    for (int v = start; v < end; v++) {
        int mx = 1;
        for (int o = 0; o < variables[v].n_outcomes; o++) {
            if (variables[v].mul_n[o] > mx) mx = variables[v].mul_n[o];
            if (variables[v].mul_d[o] > mx) mx = variables[v].mul_d[o];
        }
        mpz_mul_ui(result, result, mx);
    }
}

double p898_native(void) {
    int ps[51];
    for (int i = 0; i < 51; i++) ps[i] = 25 + i;

    build_variables(ps, 51);
    if (n_vars == 0) return 0.5;

    int mid = n_vars / 2;
    n_left = mid;

    mpz_t l_bound, r_bound, bound;
    mpz_init(l_bound); mpz_init(r_bound); mpz_init(bound);
    compute_half_bound(l_bound, 0, mid);
    compute_half_bound(r_bound, mid, n_vars);
    if (mpz_cmp(l_bound, r_bound) > 0) mpz_set(bound, l_bound);
    else mpz_set(bound, r_bound);
    shift_g = 2 * (int)mpz_sizeinbase(bound, 2) + 4;
    mpz_clear(l_bound); mpz_clear(r_bound); mpz_clear(bound);

    enumerate_right(&variables[mid], n_vars - mid, shift_g);

    for (int i = 0; i < 14; i++) { mpz_init(lnum[i]); mpz_init(lden[i]); }
    pHA = 0; cH = 0; pTA = 0; cT = 0;
    enum_left(0, 1.0, 1.0);
    for (int i = 0; i < 14; i++) { mpz_clear(lnum[i]); mpz_clear(lden[i]); }

    double tv = pHA - pTA;
    if (tv < -1.0) tv = -1.0;
    if (tv > 1.0) tv = 1.0;

    return 0.5 * (1.0 + tv);
}
