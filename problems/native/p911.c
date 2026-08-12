// Project Euler 911: Khinchin Exceptions
// Geometric mean of k_infty(rho_n) for 0 <= n <= 50.
// Uses Shallit's continued fraction recurrences (Theorems 1 and 11).
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Dynamic long long array */
typedef struct { long long *data; int len, cap; } LA;

static void la_init(LA *a) { a->cap = 16; a->len = 0; a->data = malloc(16 * sizeof(long long)); }
static void la_free(LA *a) { free(a->data); a->data = NULL; }
static void la_push(LA *a, long long v) {
    if (a->len >= a->cap) { a->cap *= 2; a->data = realloc(a->data, a->cap * sizeof(long long)); }
    a->data[a->len++] = v;
}
static long long la_last(const LA *a) { return a->data[a->len - 1]; }
static LA la_copy(const LA *a) {
    LA r; r.cap = a->cap; r.len = a->len; r.data = malloc(r.cap * sizeof(long long));
    memcpy(r.data, a->data, r.len * sizeof(long long)); return r;
}

/* Seed for rho_n continued fraction prefix.
   mode 0 = Theorem 1 (n==0), mode 1 = Theorem 11 (n>=1). */
typedef struct { LA coeffs; int mode; long long t2; } Seed;

static Seed seed_rho(int n) {
    Seed s; la_init(&s.coeffs);
    if (n == 0) {
        la_push(&s.coeffs, 0); la_push(&s.coeffs, 1); la_push(&s.coeffs, 3);
        s.mode = 0; s.t2 = 0;
        return s;
    }
    int vprime = 0;
    while ((1LL << vprime) <= n) vprime++;
    long long d = (1LL << vprime) - n;
    long long ut = 1LL << n;
    long long ud = 1LL << d;
    long long c = 0;
    for (int k = 0; k < vprime; k++)
        c += 1LL << (n - (1 << k));
    la_push(&s.coeffs, c);
    la_push(&s.coeffs, ud - 1);
    la_push(&s.coeffs, 1);
    la_push(&s.coeffs, ut - 1);
    la_push(&s.coeffs, ud);
    s.mode = 1; s.t2 = ut - 1;
    return s;
}

/* One application of Shallit Theorem 11(B) for u=2 */
static LA extend_theorem11(const LA *coeffs, long long t2) {
    LA result = la_copy(coeffs);
    long long last = la_last(coeffs);
    if (last > 1) {
        la_push(&result, t2);
        la_push(&result, 1);
        la_push(&result, last - 1);
        for (int i = coeffs->len - 2; i >= 1; i--)
            la_push(&result, coeffs->data[i]);
    } else {
        long long prelast = coeffs->data[coeffs->len - 2];
        la_push(&result, t2);
        la_push(&result, 1 + prelast);
        for (int i = coeffs->len - 3; i >= 1; i--)
            la_push(&result, coeffs->data[i]);
    }
    return result;
}

/* One application of Shallit Theorem 1(B) for u=2 */
static LA extend_theorem1(const LA *coeffs) {
    long long last = la_last(coeffs);
    LA result;
    result.cap = coeffs->cap + coeffs->len; result.len = 0;
    result.data = malloc(result.cap * sizeof(long long));
    for (int i = 0; i < coeffs->len - 1; i++)
        la_push(&result, coeffs->data[i]);
    if (last > 1) {
        la_push(&result, last + 1);
        la_push(&result, last - 1);
        for (int i = coeffs->len - 2; i >= 1; i--)
            la_push(&result, coeffs->data[i]);
    } else {
        long long prelast = coeffs->data[coeffs->len - 2];
        la_push(&result, 2 + prelast);
        for (int i = coeffs->len - 3; i >= 1; i--)
            la_push(&result, coeffs->data[i]);
    }
    return result;
}

static void canonicalize_tail(LA *cf) {
    if (cf->len > 1 && cf->data[cf->len - 1] == 1) {
        cf->data[cf->len - 2] += 1;
        cf->len--;
    }
}

static double avg_log_positive(const long long *arr, int len) {
    double s = 0.0;
    for (int i = 0; i < len; i++)
        s += log((double)arr[i]);
    return s / len;
}

static double log_khinchin_rho0(void) {
    Seed s = seed_rho(0);
    LA coeffs = s.coeffs;
    canonicalize_tail(&coeffs);

    double mu1 = 0, mu2 = 0;
    int L1 = 0, L2 = 0;

    for (int step = 1; step <= 22; step++) {
        LA next = extend_theorem1(&coeffs);
        la_free(&coeffs);
        coeffs = next;
        canonicalize_tail(&coeffs);

        if (step == 20) {
            L1 = coeffs.len - 1;
            mu1 = avg_log_positive(coeffs.data + 1, L1);
        }
        if (step == 22) {
            L2 = coeffs.len - 1;
            mu2 = avg_log_positive(coeffs.data + 1, L2);
        }
    }
    la_free(&coeffs);

    return (mu2 * L2 - mu1 * L1) / (L2 - L1);
}

static double khinchin_log_limit(int n, int max_steps) {
    if (n == 0) return log_khinchin_rho0();

    Seed s = seed_rho(n);
    LA coeffs = s.coeffs;
    long long t2 = s.t2;

    /* tail = coeffs[1:] */
    double L = (double)(coeffs.len - 1);
    double mu = 0.0;
    for (int i = 1; i < coeffs.len; i++)
        mu += log((double)coeffs.data[i]);
    mu /= L;

    long long first = coeffs.data[1];
    long long second = coeffs.data[2];
    long long last = coeffs.data[coeffs.len - 1];
    long long prelast = coeffs.data[coeffs.len - 2];

    for (int step = 0; step < max_steps; step++) {
        double mu_old = mu;
        double L_new;
        double delta;
        if (last > 1) {
            L_new = 2.0 * L + 2.0;
            delta = -log((double)last) + log((double)t2) + log((double)(last - 1));
            mu = mu * (2.0 * L / L_new) + delta / L_new;
        } else {
            L_new = 2.0 * L;
            delta = -log((double)prelast) + log((double)t2) + log((double)(prelast + 1));
            mu = mu + delta / L_new;
        }
        L = L_new;
        prelast = second;
        last = first;
        if (fabs(mu - mu_old) < 1e-15) break;
    }

    la_free(&coeffs);
    return mu;
}

double p911_native(void) {
    double total_log = 0.0;
    for (int n = 0; n <= 50; n++)
        total_log += khinchin_log_limit(n, 80);
    double avg_log = total_log / 51.0;
    return exp(avg_log);
}
