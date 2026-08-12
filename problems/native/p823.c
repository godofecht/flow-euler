#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* Project Euler 823 - Factor Shuffle.
 *
 * S(10^4, 10^16) mod 1234567891.
 *
 * Each number is represented as its sorted list of prime factors (with
 * multiplicity). A round removes the smallest prime factor from each list,
 * concatenates those removed primes into a new sorted list, appends it, and
 * drops any emptied lists.
 *
 * For large m the process reaches a regime where the k-th prime factor of the
 * newly created number becomes periodic with period k. After detecting that
 * regime we jump to extremely large m by indexing into these short cycles.
 */

#define MOD 1234567891LL

typedef struct {
    int *factors;
    int len;
    int pos;
} Pile;

static int *sieve_spf(int n) {
    int *spf = malloc((size_t)(n + 1) * sizeof(int));
    for (int i = 0; i <= n; i++) spf[i] = i;
    int limit = (int)sqrt((double)n);
    for (int i = 2; i <= limit; i++) {
        if (spf[i] == i) {
            for (long j = (long)i * i; j <= n; j += i) {
                if (spf[j] == j) spf[j] = i;
            }
        }
    }
    return spf;
}

static int count_factors(int x, int *spf) {
    int cnt = 0;
    while (x > 1) { cnt++; x /= spf[x]; }
    return cnt;
}

static void fill_factors(int x, int *spf, int *out) {
    int cnt = 0;
    while (x > 1) { out[cnt++] = spf[x]; x /= spf[x]; }
}

long long p823_native(void) {
    int n = 10000;
    long long m = 10000000000000000LL; /* 10^16 */
    long long mod = MOD;

    int *spf = sieve_spf(n);

    /* Initial piles: one per number 2..n */
    Pile *piles = malloc((size_t)n * sizeof(Pile));
    int npiles = 0;
    int total_factors = 0;

    for (int i = 2; i <= n; i++) {
        int cnt = count_factors(i, spf);
        int *f = malloc((size_t)cnt * sizeof(int));
        fill_factors(i, spf, f);
        piles[npiles].factors = f;
        piles[npiles].len = cnt;
        piles[npiles].pos = 0;
        npiles++;
        total_factors += cnt;
    }

    int k_extra = 10;
    int k_lim = (int)sqrt(2.0 * total_factors) + k_extra;

    /* Circular buffers for periodicity detection.
     * bufs[k] holds a deque of max length k. */
    int **bufs = calloc((size_t)(k_lim + 1), sizeof(int *));
    int *buf_start = calloc((size_t)(k_lim + 1), sizeof(int));
    int *buf_count = calloc((size_t)(k_lim + 1), sizeof(int));
    for (int k = 1; k <= k_lim; k++) {
        bufs[k] = malloc((size_t)k * sizeof(int));
    }

    int *count_arr = calloc((size_t)(n + 1), sizeof(int));

    int streak_needed = 2000;
    int max_rounds = 200000;
    int stable_streak = 0;
    long t = 0;
    long end_t = 0;
    int kmax = 0;
    int **patterns = NULL;

    while (t < max_rounds) {
        t++;
        int k = npiles;

        Pile *new_piles = malloc((size_t)(k + 1) * sizeof(Pile));
        int new_npiles = 0;

        /* Extract smallest remaining factor from each pile.
         * Use counting sort since all values are primes in [2, n]. */
        int max_v = 0;
        for (int idx = 0; idx < k; idx++) {
            int v = piles[idx].factors[piles[idx].pos];
            count_arr[v]++;
            if (v > max_v) max_v = v;
            piles[idx].pos++;
            if (piles[idx].pos < piles[idx].len) {
                new_piles[new_npiles++] = piles[idx];
            } else {
                free(piles[idx].factors);
            }
        }

        int *sorted_extracted = malloc((size_t)k * sizeof(int));
        int sidx = 0;
        for (int v = 2; v <= max_v; v++) {
            while (count_arr[v] > 0) {
                sorted_extracted[sidx++] = v;
                count_arr[v]--;
            }
        }

        new_piles[new_npiles].factors = sorted_extracted;
        new_piles[new_npiles].len = k;
        new_piles[new_npiles].pos = 0;
        new_npiles++;

        free(piles);
        piles = new_piles;
        npiles = new_npiles;

        /* Update per-column buffers and check periodicity. */
        if (t <= k_lim) {
            int mm = k < k_lim ? k : k_lim;
            for (int kk = 1; kk <= mm; kk++) {
                int v = sorted_extracted[kk - 1];
                int idx = (buf_start[kk] + buf_count[kk]) % kk;
                if (buf_count[kk] == kk) buf_start[kk] = (buf_start[kk] + 1) % kk;
                else buf_count[kk]++;
                bufs[kk][idx] = v;
            }
            for (int kk = mm + 1; kk <= k_lim; kk++) {
                int idx = (buf_start[kk] + buf_count[kk]) % kk;
                if (buf_count[kk] == kk) buf_start[kk] = (buf_start[kk] + 1) % kk;
                else buf_count[kk]++;
                bufs[kk][idx] = 1;
            }
            stable_streak = 0;
            continue;
        }

        int all_ok = 1;
        int mm = k < k_lim ? k : k_lim;

        for (int kk = 1; kk <= mm; kk++) {
            int v = sorted_extracted[kk - 1];
            if (bufs[kk][buf_start[kk]] != v) all_ok = 0;
            int idx = (buf_start[kk] + buf_count[kk]) % kk;
            if (buf_count[kk] == kk) buf_start[kk] = (buf_start[kk] + 1) % kk;
            else buf_count[kk]++;
            bufs[kk][idx] = v;
        }
        for (int kk = mm + 1; kk <= k_lim; kk++) {
            if (bufs[kk][buf_start[kk]] != 1) all_ok = 0;
            int idx = (buf_start[kk] + buf_count[kk]) % kk;
            if (buf_count[kk] == kk) buf_start[kk] = (buf_start[kk] + 1) % kk;
            else buf_count[kk]++;
            bufs[kk][idx] = 1;
        }

        if (all_ok) {
            stable_streak++;
            if (stable_streak >= streak_needed) {
                end_t = t;
                patterns = calloc((size_t)(k_lim + 1), sizeof(int *));
                kmax = 0;
                for (int kk = 1; kk <= k_lim; kk++) {
                    patterns[kk] = malloc((size_t)kk * sizeof(int));
                    for (int j = 0; j < kk; j++) {
                        patterns[kk][j] = bufs[kk][(buf_start[kk] + j) % kk];
                    }
                    for (int j = 0; j < kk; j++) {
                        if (patterns[kk][j] != 1) { kmax = kk; break; }
                    }
                }
                break;
            }
        } else {
            stable_streak = 0;
        }
    }

    if (patterns == NULL) {
        fprintf(stderr, "p823: periodicity not detected within %d rounds\n", max_rounds);
        return -1;
    }

    /* In the periodic regime, x(u, k) = patterns[k][(u - end_t - 1) % k].
     * After m rounds, the number added at round (m - d) has been divided d
     * times, leaving the suffix product of its factors from index (d+1) on. */
    __int128 r0 = (__int128)m - (__int128)end_t - 1;
    long long total = 0;

    for (int d = 0; d < kmax; d++) {
        __int128 r = r0 - d;
        int dk = d + 1;
        if (patterns[dk][(long long)(r % dk)] == 1) continue;

        long long prod = 1;
        for (int kk = kmax; kk > d; kk--) {
            int v = patterns[kk][(long long)(r % kk)];
            if (v != 1) prod = (prod * v) % mod;
        }
        total = (total + prod) % mod;
    }

    return total;
}
