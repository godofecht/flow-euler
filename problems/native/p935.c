/* Project Euler 935: Rolling Square
 *
 * Count b in (0,1) such that the rolling square returns to its initial
 * position within N rolls. F(10^8).
 *
 * Uses Mobius inversion with the Dirichlet hyperbola method,
 * accelerated by a Du Jiao sieve for the Mertens function.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef int64_t i64;

/* ---- Linear sieve for Mobius ---- */

static int8_t *mu_arr;
static i64 *M_small;     /* Mertens function prefix sums */
static i64 *Modd_small;  /* odd-restricted Mertens prefix sums */
static i64 sieve_limit;

static void build_mu_prefix(i64 limit) {
    sieve_limit = limit;
    mu_arr = (int8_t *)calloc(limit + 1, sizeof(int8_t));
    int8_t *is_comp = (int8_t *)calloc(limit + 1, sizeof(int8_t));
    i64 *primes = (i64 *)malloc(limit * sizeof(i64) / 4 + 100);
    i64 nprimes = 0;

    mu_arr[1] = 1;
    for (i64 i = 2; i <= limit; i++) {
        if (!is_comp[i]) {
            primes[nprimes++] = i;
            mu_arr[i] = -1;
        }
        for (i64 j = 0; j < nprimes; j++) {
            i64 ip = i * primes[j];
            if (ip > limit) break;
            is_comp[ip] = 1;
            if (i % primes[j] == 0) {
                mu_arr[ip] = 0;
                break;
            }
            mu_arr[ip] = -mu_arr[i];
        }
    }

    M_small = (i64 *)calloc(limit + 1, sizeof(i64));
    Modd_small = (i64 *)calloc(limit + 1, sizeof(i64));

    i64 s = 0, so = 0;
    for (i64 i = 1; i <= limit; i++) {
        s += mu_arr[i];
        M_small[i] = s;
        if (i & 1) so += mu_arr[i];
        Modd_small[i] = so;
    }

    free(is_comp);
    free(primes);
}

/* ---- Du Jiao sieve for Mertens function ---- */

/* Hash table for memoized M(n) and Modd(n) */
typedef struct { i64 n; i64 val; } MemoEntry;

#define MEMO_SIZE (1 << 20)
#define MEMO_MASK (MEMO_SIZE - 1)

static MemoEntry *memo_M;
static MemoEntry *memo_Modd;

static i64 memo_lookup(MemoEntry *table, i64 n) {
    i64 h = (i64)((uint64_t)n * 2654435761ULL) & MEMO_MASK;
    for (int i = 0; i < MEMO_SIZE; i++) {
        i64 idx = (h + i) & MEMO_MASK;
        if (table[idx].n == n && table[idx].val != 0) return table[idx].val;
        if (table[idx].val == 0 && table[idx].n == 0) return INT64_MIN; /* not found */
    }
    return INT64_MIN;
}

static void memo_insert(MemoEntry *table, i64 n, i64 val) {
    i64 h = (i64)((uint64_t)n * 2654435761ULL) & MEMO_MASK;
    for (int i = 0; i < MEMO_SIZE; i++) {
        i64 idx = (h + i) & MEMO_MASK;
        if (table[idx].val == 0) {
            table[idx].n = n;
            table[idx].val = val;
            return;
        }
    }
}

static i64 M_func(i64 n);
static i64 Modd_func(i64 n);

static i64 M_func(i64 n) {
    if (n <= sieve_limit) return M_small[n];
    i64 found = memo_lookup(memo_M, n);
    if (found != INT64_MIN) return found;

    i64 res = 1;
    i64 l = 2;
    while (l <= n) {
        i64 q = n / l;
        i64 r = n / q;
        res -= (r - l + 1) * M_func(q);
        l = r + 1;
    }

    memo_insert(memo_M, n, res);
    return res;
}

static i64 Modd_func(i64 n) {
    if (n <= 0) return 0;
    if (n <= sieve_limit) return Modd_small[n];
    i64 found = memo_lookup(memo_Modd, n);
    if (found != INT64_MIN) return found;

    i64 res = M_func(n) + Modd_func(n / 2);
    memo_insert(memo_Modd, n, res);
    return res;
}

/* ---- Helper functions ---- */

static i64 tri(i64 n) { return n * (n + 1) / 2; }

static i64 sum_mu_odd(i64 l, i64 r) {
    return Modd_func(r) - Modd_func(l - 1);
}

static i64 sum_mu_2mod4(i64 l, i64 r) {
    return -Modd_func(r / 2) + Modd_func((l - 1) / 2);
}

/* ---- Class sum ---- */

static i64 class_sum(i64 X, int cls) {
    i64 A = 0, B = 0;
    i64 l = 1;
    while (l <= X) {
        i64 q = X / l;
        i64 r = X / q;

        i64 odd_mu = sum_mu_odd(l, r);
        i64 mu2 = sum_mu_2mod4(l, r);

        if (cls == 0) { /* div4 */
            A += q * ((q / 4) * odd_mu + (q / 2) * mu2);
            B += (4 * tri(q / 4)) * odd_mu + (2 * tri(q / 2)) * mu2;
        } else if (cls == 1) { /* 2mod4 */
            i64 c = (q + 2) / 4;
            A += q * (c * odd_mu + ((q + 1) / 2) * mu2);
            B += (2 * c * c) * odd_mu + (tri(q) - 2 * tri(q / 2)) * mu2;
        } else { /* odd */
            A += q * (((q + 1) / 2) * odd_mu);
            B += (tri(q) - 2 * tri(q / 2)) * odd_mu;
        }

        l = r + 1;
    }
    return A - B;
}

/* ---- F(N) ---- */

static i64 F_935(i64 N, i64 max_N) {
    if (N < 0) return 0;

    i64 X1 = N + 1;
    i64 X2 = N / 2 + 1;
    i64 X3 = N / 4 + 1;

    i64 res = 0;
    res += class_sum(X1, 0);  /* div4 */
    res += class_sum(X2, 1);  /* 2mod4 */
    res += class_sum(X3, 2);  /* odd */

    /* Remove u=1 (odd class) */
    res -= X3 - 1;

    /* Integer L = h cases: 4*(h-1) rolls, h >= 2 */
    res += N / 4;

    return res;
}

long long p935_native(void) {
    i64 N = 1;
    for (int i = 0; i < 8; i++) N *= 10;

    i64 max_X = N + 1;
    i64 limit = (i64)pow((double)max_X, 2.0 / 3.0) + 10;
    build_mu_prefix(limit);

    memo_M = (MemoEntry *)calloc(MEMO_SIZE, sizeof(MemoEntry));
    memo_Modd = (MemoEntry *)calloc(MEMO_SIZE, sizeof(MemoEntry));

    i64 result = F_935(N, N);

    free(mu_arr); free(M_small); free(Modd_small);
    free(memo_M); free(memo_Modd);

    return result;
}
