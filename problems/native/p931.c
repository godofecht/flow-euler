#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Project Euler 931: Totient Graph
   Compute T(10^12) mod 715827883.
   Uses Min25-style prime counting/sum sieve and a closed-form summation.
   Port of the Python reference solver. */

typedef int64_t i64;
typedef __int128 i128;

#define MOD 715827883LL
#define TARGET 1000000000000LL

/* ---------- Sieve primes (odd-only) ---------- */
static int *sieve_primes(int limit, int *count) {
    if (limit < 2) { *count = 0; return NULL; }
    if (limit == 2) { int *p = malloc(sizeof(int)); p[0] = 2; *count = 1; return p; }

    int size = limit / 2 + 1;
    char *is_comp = calloc(size, 1);
    int *primes = malloc((size + 1) * sizeof(int));
    int pc = 0;
    primes[pc++] = 2;

    int r = (int)sqrt((double)limit);
    for (int i = 1; i <= r / 2; i++) {
        if (!is_comp[i]) {
            int p = 2 * i + 1;
            int start = (p * p) / 2;
            int step = p;
            for (int j = start; j < size; j += step)
                is_comp[j] = 1;
        }
    }
    for (int i = 1; i < size; i++)
        if (!is_comp[i]) primes[pc++] = 2 * i + 1;

    free(is_comp);
    *count = pc;
    return primes;
}

/* ---------- F(x) = x(x+1)/2 mod mod ---------- */
static i64 tri(i64 x, i64 mod) {
    /* x up to 10^12, x*(x+1) up to 10^24, need i128 */
    i128 val = (i128)x * (x + 1) / 2;
    return (i64)(val % mod);
}

/* ---------- Min25 prime count and sum sieve ---------- */
static int v;
static int *primes_arr;
static int num_primes;
static i64 *vals;     /* vals[i] = n//i for 1<=i<=v */
static i64 *g_small;  /* pi(x) for x <= v */
static i64 *h_small;  /* prime sum mod for x <= v */
static i64 *g_large;  /* pi(n//i) for 1<=i<=v */
static i64 *h_large;  /* prime sum mod for n//i */

static void min25_init(i64 n, i64 mod) {
    v = (int)sqrt((double)n);
    primes_arr = sieve_primes(v, &num_primes);

    vals = calloc(v + 1, sizeof(i64));
    g_small = calloc(v + 1, sizeof(i64));
    h_small = calloc(v + 1, sizeof(i64));
    g_large = calloc(v + 1, sizeof(i64));
    h_large = calloc(v + 1, sizeof(i64));

    /* Initialize small */
    for (int x = 0; x <= v; x++) {
        if (x >= 2) {
            g_small[x] = x - 1;
            h_small[x] = ((i128)x * (x + 1) / 2 - 1) % mod;
        } else {
            g_small[x] = 0;
            h_small[x] = 0;
        }
    }

    /* Initialize large */
    for (int i = 1; i <= v; i++) {
        i64 x = n / i;
        vals[i] = x;
        if (x >= 2) {
            g_large[i] = x - 1;
            h_large[i] = ((i128)x * (x + 1) / 2 - 1) % mod;
        } else {
            g_large[i] = 0;
            h_large[i] = 0;
        }
    }

    /* Process primes */
    for (int pi = 0; pi < num_primes; pi++) {
        i64 p = primes_arr[pi];
        i64 p2 = p * p;
        if (p2 > n) break;

        i64 g_p1 = g_small[p - 1];
        i64 h_p1 = h_small[p - 1];

        /* Update large values, increasing i */
        i64 i_max = n / p2;
        if (i_max > v) i_max = v;

        for (int i = 1; i <= i_max; i++) {
            i64 y = vals[i] / p;
            i64 g_y, h_y;
            if (y <= v) {
                g_y = g_small[y];
                h_y = h_small[y];
            } else {
                int k = (int)(n / y);
                g_y = g_large[k];
                h_y = h_large[k];
            }

            g_large[i] -= g_y - g_p1;

            i64 diff = h_y - h_p1;
            if (diff < 0) diff += mod;
            h_large[i] = (h_large[i] - (i128)p * diff % mod) % mod;
            if (h_large[i] < 0) h_large[i] += mod;
        }

        /* Update small values descending */
        for (int x = v; x >= p2; x--) {
            int y = x / (int)p;
            g_small[x] -= g_small[y] - g_p1;

            i64 diff = h_small[y] - h_p1;
            if (diff < 0) diff += mod;
            h_small[x] = (h_small[x] - (i128)p * diff % mod) % mod;
            if (h_small[x] < 0) h_small[x] += mod;
        }
    }
}

static i64 pi(i64 x) {
    if (x <= v) return g_small[x];
    return g_large[TARGET / x];
}

static i64 psum(i64 x) {
    if (x <= v) return h_small[x];
    return h_large[TARGET / x];
}

static i64 compute_T_mod(i64 N, i64 mod) {
    if (N <= 1) return 0;

    int sqrtN = (int)sqrt((double)N);
    min25_init(N, mod);

    i64 total = 0;

    /* Part 1: primes p <= sqrt(N), all exponents e >= 1 */
    for (int pi = 0; pi < num_primes; pi++) {
        i64 p = primes_arr[pi];
        if (p * p > N) break;

        /* e = 1 */
        i64 x1 = N / p;
        i64 x2 = N / (p * p);
        i64 f = (tri(x1, mod) - (p % mod) * tri(x2, mod)) % mod;
        if (f < 0) f += mod;
        total = (total + ((p - 2) % mod) * f) % mod;

        /* e >= 2 */
        i64 pe = p * p;
        i64 p_pow = p; /* p^{e-1} for current e, starts at e=2 */
        while (pe <= N) {
            i64 xe = N / pe;
            i64 xnext = N / (pe * p);
            i64 A = (p - 1) * p_pow - 1;
            i64 f2 = (tri(xe, mod) - (p % mod) * tri(xnext, mod)) % mod;
            if (f2 < 0) f2 += mod;
            total = (total + (A % mod) * f2) % mod;
            p_pow *= p;
            pe *= p;
        }
    }

    /* Part 2: primes p > sqrt(N), only e=1. Group by q=floor(N/p) */
    for (i64 q = 1; q < sqrtN; q++) {
        i64 hi = N / q;
        if (hi <= sqrtN) break;

        i64 lo = N / (q + 1) + 1;
        if (lo <= sqrtN) lo = sqrtN + 1;
        if (lo > hi) continue;

        i64 sum_p = (psum(hi) - psum(lo - 1)) % mod;
        if (sum_p < 0) sum_p += mod;
        i64 cnt_p = pi(hi) - pi(lo - 1);
        i64 term2 = (sum_p - (2 * (cnt_p % mod)) % mod) % mod;
        if (term2 < 0) term2 += mod;
        total = (total + tri(q, mod) * term2) % mod;
    }

    total %= mod;
    if (total < 0) total += mod;

    /* Cleanup */
    free(vals); free(g_small); free(h_small);
    free(g_large); free(h_large);
    free(primes_arr);

    return total;
}

long long p931_native(void) {
    return compute_T_mod(TARGET, MOD);
}
