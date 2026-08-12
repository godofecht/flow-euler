/* Project Euler 962: Counting Triangles
 *
 * count_triangles(10^6) = 7259046
 *
 * Counts triangles (a,b,c) with a<=b<=c<a+b and a+b+c<=limit such that
 * a^3*(a+b-c)*(a+b+c) / (b*(a+b)^2) is a perfect square.
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef long long i64;

static i64 isqrt_i64(i64 n) {
    if (n <= 0) return 0;
    i64 r = (i64)sqrt((double)n);
    while (r > 0 && r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r;
}

static i64 isqrt_ceil_i64(i64 n) {
    if (n <= 0) return 0;
    i64 r = isqrt_i64(n);
    return (r * r == n) ? r : r + 1;
}

static i64 icbrt_floor(i64 n) {
    if (n <= 0) return 0;
    i64 x = (i64)cbrt((double)n);
    while ((x + 1) * (x + 1) * (x + 1) <= n) x++;
    while (x * x * x > n) x--;
    return x;
}

static i64 igcd(i64 a, i64 b) {
    while (b) { i64 t = a % b; a = b; b = t; }
    return a;
}

static i64 count_with_parity(i64 lo, i64 hi, i64 parity) {
    if (hi < lo) return 0;
    if ((lo & 1) != parity) lo++;
    if (lo > hi) return 0;
    return (hi - lo) / 2 + 1;
}

static int *spf_sieve(int limit) {
    int *spf = malloc((size_t)(limit + 1) * sizeof(int));
    for (int i = 0; i <= limit; i++) spf[i] = i;
    int sq = (int)sqrt((double)limit);
    for (int i = 2; i <= sq; i++) {
        if (spf[i] == i) {
            for (int j = i * i; j <= limit; j += i) {
                if (spf[j] == j) spf[j] = i;
            }
        }
    }
    return spf;
}

static int *compute_squarefree(int limit, int *spf) {
    int *sf = malloc((size_t)(limit + 1) * sizeof(int));
    sf[0] = 0; sf[1] = 1;
    for (int n = 2; n <= limit; n++) {
        int p = spf[n];
        int m = n / p;
        if (m % p == 0)
            sf[n] = sf[m / p];
        else
            sf[n] = p * sf[m];
    }
    return sf;
}

static i64 count_triangles(i64 limit) {
    i64 k_limit = icbrt_floor(2 * limit * limit) + 2;
    i64 factor_limit = (limit / 3 + 1 > k_limit + 1)
                       ? limit / 3 + 1 : k_limit + 1;

    int *spf = spf_sieve((int)factor_limit);
    int *squarefree = compute_squarefree((int)factor_limit, spf);

    i64 total = 0;
    for (i64 k = 2; k <= k_limit; k++) {
        i64 max_s = limit / k;
        if (max_s == 0) break;
        i64 max_product = max_s * max_s;
        for (i64 v = (k + 1) / 2; v < k; v++) {
            if (igcd(v, k) != 1) continue;
            i64 u = k - v;
            i64 d = squarefree[(int)u];
            if (v * d > max_product) continue;

            i64 r_limit = u * max_s / (u + 2 * v);
            if (r_limit <= 0) continue;

            i64 common_dv = igcd(d, v);
            i64 d_reduced = d / common_dv;
            i64 v_reduced = v / common_dv;
            i64 numerator = u + 2 * v;

            for (i64 r = 1; r <= r_limit; r++) {
                i64 common_rv = igcd(r, v_reduced);
                i64 remaining_r = r / common_rv;
                i64 numerator_kernel = squarefree[(int)remaining_r];
                i64 common = igcd(d_reduced, numerator_kernel);
                i64 base = (v_reduced / common_rv) *
                           (d_reduced * numerator_kernel / (common * common));
                if (base > max_s) continue;

                i64 s_min = (numerator * r + u - 1) / u;
                i64 n_min = isqrt_ceil_i64((s_min + base - 1) / base);
                i64 n_max = isqrt_i64(max_s / base);
                if (base & 1) {
                    total += count_with_parity(n_min, n_max, r & 1);
                } else if ((r & 1) == 0) {
                    i64 cnt = n_max - n_min + 1;
                    if (cnt > 0) total += cnt;
                }
            }
        }
    }

    free(spf);
    free(squarefree);
    return total;
}

long long p962_native(void) {
    return count_triangles(1000000);
}
