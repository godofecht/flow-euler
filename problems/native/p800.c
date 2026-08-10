#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef long long i64;

static int is_prime_small(i64 x) {
    if (x <= 1) return 0;
    if (x <= 3) return 1;
    if ((x & 1) == 0) return 0;
    for (i64 i = 3; i * i <= x; i += 2)
        if (x % i == 0) return 0;
    return 1;
}

long long pe_solve(void) {
    i64 a = 800800, b = 800800;
    double limit = log10((double)a) * (double)b;
    i64 maxn = (i64)(limit / log10(2.0)) + 10000;
    unsigned char *comp = calloc((size_t)maxn + 1, 1);
    for (i64 i = 2; i * i <= maxn; i++)
        if (!comp[i])
            for (i64 j = i * i; j <= maxn; j += i) comp[j] = 1;
    int *primes = malloc((size_t)(maxn / 10) * sizeof(int));
    int pc = 0;
    for (i64 i = 2; i <= maxn; i++)
        if (!comp[i]) primes[pc++] = (int)i;
    int *pp = malloc((size_t)(maxn + 2) * sizeof(int));
    int p_index = 0;
    for (i64 x = 1; x <= maxn; x++) {
        while (p_index < pc && primes[p_index] <= x) p_index++;
        pp[x] = p_index; /* number of primes <= x */
    }
    i64 count = 0;
    /* iterate primes p */
    for (int pi = 0; pi < pc; pi++) {
        i64 p = primes[pi];
        double p_log = log10((double)p);
        if (2.0 * (double)p * p_log > limit) break;
        i64 lo = (pi + 1 < pc) ? primes[pi + 1] : p + 1;
        while (lo <= maxn && comp[lo]) lo++;
        i64 hi = (i64)(limit / p_log);
        if (hi > maxn) hi = maxn;
        while (lo < hi) {
            i64 q = (lo + hi) / 2;
            double q_log = log10((double)q);
            if ((double)p * q_log + (double)q * p_log > limit)
                hi = q;
            else
                lo = q + 1;
        }
        i64 q = lo;
        if (q > maxn || comp[q]) {
            while (q <= maxn && comp[q]) q++;
            if (q > maxn) break;
        }
        double q_log = log10((double)q);
        double val = (double)p * q_log + (double)q * p_log;
        if (fabs(val - limit) < 1e-9)
            count += pp[q] - pp[p];
        else if (val > limit)
            count += pp[q] - pp[p] - 1;
        else {
            q++;
            while (q <= maxn && comp[q]) q++;
            if (q > maxn) break;
            count += pp[q] - pp[p] - 1;
        }
    }
    free(comp);
    free(primes);
    free(pp);
    (void)is_prime_small;
    return count;
}
