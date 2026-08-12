/* Project Euler 893: Matchsticks
 *
 * T(10^6) = sum_{n=1..10^6} M(n), where M(n) is the minimum matchsticks
 * to display an expression (digits, +, x with standard precedence)
 * evaluating to n.
 *
 * Phases:
 *  1. P(n): min cost using only digits and multiplication.
 *  2. M(n): min cost allowing two-term sums a+b.
 *  3. Three-term fix: allow a+b+c when it beats 1- and 2-term forms.
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef long long i64;

static const int DIGIT_COST[10] = {6, 2, 5, 5, 4, 5, 6, 3, 7, 6};
#define OP_COST 2
#define MAXCOST 48

static int *bucket_data;
static int bucket_off[MAXCOST + 1];

static int cmp_int(const void *a, const void *b) {
    return *(const int *)a - *(const int *)b;
}

static int *spf_sieve(int n_max) {
    int *spf = malloc((size_t)(n_max + 1) * sizeof(int));
    for (int i = 0; i <= n_max; i++) spf[i] = i;
    int limit = (int)sqrt((double)n_max);
    for (int i = 2; i <= limit; i++) {
        if (spf[i] == i) {
            for (int j = i * i; j <= n_max; j += i) {
                if (spf[j] == j) spf[j] = i;
            }
        }
    }
    return spf;
}

static unsigned char *compute_P(int n_max, int *spf) {
    unsigned char *P = malloc((size_t)(n_max + 1));
    P[0] = 0;
    for (int n = 1; n <= n_max; n++)
        P[n] = P[n / 10] + DIGIT_COST[n % 10];
    P[0] = 255;

    int divs[4096];
    int fp[16], fe[16];

    for (int n = 2; n <= n_max; n++) {
        if (spf[n] == n) continue;

        int m = n, nf = 0;
        while (m > 1) {
            int p = spf[m], e = 0;
            while (m % p == 0) { m /= p; e++; }
            fp[nf] = p; fe[nf] = e; nf++;
        }

        int ndiv = 1;
        divs[0] = 1;
        for (int i = 0; i < nf; i++) {
            int p = fp[i], e = fe[i], prev = ndiv, pe = 1;
            for (int j = 0; j < e; j++) {
                pe *= p;
                for (int k = 0; k < prev; k++)
                    divs[ndiv++] = divs[k] * pe;
            }
        }

        int r = (int)sqrt((double)n);
        int best = P[n];
        for (int i = 0; i < ndiv; i++) {
            int d = divs[i];
            if (d > 1 && d <= r) {
                int q = n / d;
                int cand = P[d] + P[q] + OP_COST;
                if (cand < best) best = cand;
            }
        }
        P[n] = best;
    }
    return P;
}

static void build_buckets(int n_max, unsigned char *P) {
    int counts[MAXCOST + 1] = {0};
    for (int x = 1; x <= n_max; x++) {
        int c = P[x];
        if (c <= MAXCOST) counts[c]++;
    }
    bucket_off[0] = 0;
    for (int c = 1; c <= MAXCOST; c++)
        bucket_off[c] = bucket_off[c - 1] + counts[c - 1];
    bucket_data = malloc((size_t)n_max * sizeof(int));
    int pos[MAXCOST + 1];
    memcpy(pos, bucket_off, sizeof(pos));
    for (int x = 1; x <= n_max; x++) {
        int c = P[x];
        if (c <= MAXCOST) bucket_data[pos[c]++] = x;
    }
}

/* Enumerate pairs from buckets[c1] x buckets[c2] with a+b <= n_max.
 * For c1==c2 only pairs with i<=j.  Updates target[s] = min(target[s], cost). */
static void enum_pairs(int n_max, int c1, int c2, int cost,
                        unsigned char *target) {
    if (c1 > MAXCOST || c2 > MAXCOST) return;
    int *A = bucket_data + bucket_off[c1];
    int alen = bucket_off[c1 + 1] - bucket_off[c1];
    int *B = bucket_data + bucket_off[c2];
    int blen = bucket_off[c2 + 1] - bucket_off[c2];
    if (alen == 0 || blen == 0) return;

    if (c1 == c2) {
        for (int i = 0; i < alen; i++) {
            int a = A[i];
            int max_b = n_max - a;
            int hi = alen - 1;
            while (hi >= i && A[hi] > max_b) hi--;
            for (int j = i; j <= hi; j++) {
                int s = a + A[j];
                if (cost < target[s]) target[s] = cost;
            }
        }
    } else {
        int *outer, *inner, olen, ilen;
        if (alen <= blen) {
            outer = A; olen = alen; inner = B; ilen = blen;
        } else {
            outer = B; olen = blen; inner = A; ilen = alen;
        }
        int hi = ilen - 1;
        for (int oi = 0; oi < olen; oi++) {
            int a = outer[oi];
            int max_b = n_max - a;
            while (hi >= 0 && inner[hi] > max_b) hi--;
            for (int j = 0; j <= hi; j++) {
                int s = a + inner[j];
                if (cost < target[s]) target[s] = cost;
            }
        }
    }
}

static void compute_M_two_terms(int n_max, unsigned char *P, unsigned char *M) {
    memcpy(M, P, (size_t)(n_max + 1));

    for (int c1 = 2; c1 <= 32; c1++) {
        int max_c2 = 32 - c1;
        if (max_c2 < c1) break;
        for (int c2 = c1; c2 <= max_c2; c2++) {
            int cand_cost = c1 + c2 + OP_COST;
            enum_pairs(n_max, c1, c2, cand_cost, M);
        }
    }

    /* Fallback scan for high-cost cases using low P(a) summands. */
    int low_count = 0;
    for (int c = 2; c <= 18; c++)
        low_count += bucket_off[c + 1] - bucket_off[c];
    int *low_sum = malloc((size_t)low_count * sizeof(int));
    int li = 0;
    for (int c = 2; c <= 18; c++)
        for (int j = bucket_off[c]; j < bucket_off[c + 1]; j++)
            low_sum[li++] = bucket_data[j];
    qsort(low_sum, (size_t)low_count, sizeof(int), cmp_int);

    for (int n = 1; n <= n_max; n++) {
        if (M[n] > 34) {
            int best = M[n];
            int half = n / 2;
            for (int j = 0; j < low_count && low_sum[j] <= half; j++) {
                int a = low_sum[j];
                int b = n - a;
                int cand = P[a] + P[b] + OP_COST;
                if (cand < best) best = cand;
            }
            M[n] = best;
        }
    }
    free(low_sum);
}

static void apply_three_term_fix(int n_max, unsigned char *P, unsigned char *M) {
    unsigned char *best2 = malloc((size_t)(n_max + 1));
    memset(best2, 255, (size_t)(n_max + 1));

    for (int c1 = 2; c1 <= 30; c1++) {
        int max_c2 = 30 - c1;
        if (max_c2 < c1) break;
        for (int c2 = c1; c2 <= max_c2; c2++) {
            int psum = c1 + c2;
            enum_pairs(n_max, c1, c2, psum, best2);
        }
    }

    /* Cheap third terms: P(c) <= 20 */
    int cheap_count = 0;
    for (int c = 2; c <= 20; c++)
        cheap_count += bucket_off[c + 1] - bucket_off[c];
    int *cheap_c = malloc((size_t)cheap_count * sizeof(int));
    int ci = 0;
    for (int c = 2; c <= 20; c++)
        for (int j = bucket_off[c]; j < bucket_off[c + 1]; j++)
            cheap_c[ci++] = bucket_data[j];
    qsort(cheap_c, (size_t)cheap_count, sizeof(int), cmp_int);

    /* Small remainders r where best2[r] <= 8 */
    int small_count = 0;
    for (int r = 2; r <= n_max; r++)
        if (best2[r] != 255 && best2[r] <= 8) small_count++;
    int *small_r = malloc((size_t)small_count * sizeof(int));
    int si = 0;
    for (int r = 2; r <= n_max; r++)
        if (best2[r] != 255 && best2[r] <= 8) small_r[si++] = r;

    for (int n = 1; n <= n_max; n++) {
        int cur = M[n];
        if (cur < 33) continue;

        int target = cur - 1;
        int psum_limit = target - 4;
        if (psum_limit < 6) continue;

        int best = cur;

        for (int j = 0; j < cheap_count; j++) {
            int c = cheap_c[j];
            if (c >= n) break;
            int pc = P[c];
            if (pc > psum_limit) continue;
            int r = n - c;
            int v = best2[r];
            if (v != 255 && v + pc <= psum_limit) {
                int cand = v + pc + 4;
                if (cand < best) {
                    best = cand;
                    if (best == target) break;
                }
            }
        }

        if (best > target) {
            for (int j = 0; j < small_count; j++) {
                int r = small_r[j];
                if (r >= n) break;
                int v = best2[r];
                int c = n - r;
                int pc = P[c];
                if (v + pc <= psum_limit) {
                    int cand = v + pc + 4;
                    if (cand < best) {
                        best = cand;
                        if (best == target) break;
                    }
                }
            }
        }

        M[n] = best;
    }

    free(best2);
    free(cheap_c);
    free(small_r);
}

long long p893_native(void) {
    int N = 1000000;

    int *spf = spf_sieve(N);
    unsigned char *P = compute_P(N, spf);
    free(spf);

    build_buckets(N, P);

    unsigned char *M = malloc((size_t)(N + 1));
    compute_M_two_terms(N, P, M);
    apply_three_term_fix(N, P, M);

    i64 total = 0;
    for (int n = 1; n <= N; n++)
        total += M[n];

    free(P);
    free(M);
    free(bucket_data);

    return total;
}
