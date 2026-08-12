#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Project Euler 834 - Add and Divide.
 *
 * U(N) = sum_{n=3..N} T(n), where T(n) is the sum of indices m such that
 * the m-th term of the sequence is divisible by (n+m).
 *
 * For each n, let even = n (if n even) or n-1 (if n odd), and B = n-1 or n.
 * A = odd part of even, two = 2^v2(even).
 * Valid divisors d of A*B come in two forms: odd d|A*B with d > n, or d = two*p
 * with p|A*B and two*p > n. T(n) = sum of (d - n) over all valid d.
 *
 * Rolling cache: for consecutive n, one of A or B repeats.
 */

typedef long long i64;

static int *spf;

static void build_spf(int limit) {
    spf = (int *)calloc((size_t)(limit + 1), sizeof(int));
    int *primes = (int *)malloc((size_t)(limit + 1) * sizeof(int));
    int pcnt = 0;
    for (int i = 2; i <= limit; i++) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes[pcnt++] = i;
        }
        for (int j = 0; j < pcnt; j++) {
            int ip = i * primes[j];
            if (ip > limit) break;
            spf[ip] = primes[j];
            if (primes[j] == spf[i]) break;
        }
    }
    free(primes);
    spf[0] = 0;
    if (limit >= 1) spf[1] = 1;
}

/* Divisor data: sorted divisors, prefix sums, total sum */
typedef struct {
    int *divs;
    i64 *pre;  /* prefix sums, pre[0]=0, pre[i] = sum of divs[0..i-1] */
    int count;
    i64 total;
} DivData;

static int cmp_int(const void *a, const void *b) {
    return *(const int *)a - *(const int *)b;
}

static DivData divisors_with_prefix(int odd_x) {
    DivData dd;
    if (odd_x == 1) {
        dd.divs = (int *)malloc(sizeof(int));
        dd.divs[0] = 1;
        dd.pre = (i64 *)malloc(2 * sizeof(i64));
        dd.pre[0] = 0;
        dd.pre[1] = 1;
        dd.count = 1;
        dd.total = 1;
        return dd;
    }

    /* Factorize using SPF */
    int factors_p[64], factors_e[64], nfac = 0;
    int x = odd_x;
    while (x > 1) {
        int p = spf[x];
        int e = 0;
        while (x % p == 0) { x /= p; e++; }
        factors_p[nfac] = p;
        factors_e[nfac] = e;
        nfac++;
    }

    /* Generate divisors */
    int cap = 1;
    for (int i = 0; i < nfac; i++) cap *= (factors_e[i] + 1);
    int *divs = (int *)malloc((size_t)cap * sizeof(int));
    int ndiv = 1;
    divs[0] = 1;
    for (int i = 0; i < nfac; i++) {
        int p = factors_p[i], e = factors_e[i];
        int base_ndiv = ndiv;
        int mult = 1;
        for (int k = 1; k <= e; k++) {
            mult *= p;
            for (int j = 0; j < base_ndiv; j++) {
                divs[ndiv++] = divs[j] * mult;
            }
        }
    }

    qsort(divs, (size_t)ndiv, sizeof(int), cmp_int);

    i64 *pre = (i64 *)malloc((size_t)(ndiv + 1) * sizeof(i64));
    i64 s = 0;
    for (int i = 0; i < ndiv; i++) {
        s += divs[i];
        pre[i + 1] = s;
    }

    dd.divs = divs;
    dd.pre = pre;
    dd.count = ndiv;
    dd.total = s;
    return dd;
}

static void free_divdata(DivData *dd) {
    free(dd->divs);
    free(dd->pre);
}

/* Binary search: rightmost index where divs[idx] <= bound */
static int bisect_right(int *divs, int count, int bound) {
    int lo = 0, hi = count;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (divs[mid] <= bound) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

/* Over pairs (a in small, b in large), count and sum where a*b > bound */
static void count_sum_products_gt(
    int *small, int ns, int *large, int nl, i64 *lpre, i64 ltotal,
    int bound, i64 *cnt_out, i64 *sum_out)
{
    i64 cnt = 0, s = 0;
    for (int i = 0; i < ns; i++) {
        int a = small[i];
        int idx = bisect_right(large, nl, bound / a);
        if (idx != nl) {
            cnt += nl - idx;
            s += (i64)a * (ltotal - lpre[idx]);
        }
    }
    *cnt_out = cnt;
    *sum_out = s;
}

static i64 compute_T_from_divdata(i64 n, int two, DivData *dA, DivData *dB) {
    int *small, *large;
    int ns, nl;
    i64 *preL, sumL;

    if (dA->count <= dB->count) {
        small = dA->divs; ns = dA->count;
        large = dB->divs; nl = dB->count;
        preL = dB->pre; sumL = dB->total;
    } else {
        small = dB->divs; ns = dB->count;
        large = dA->divs; nl = dA->count;
        preL = dA->pre; sumL = dA->total;
    }

    i64 cnt1, sum1, cnt2, sum2;
    count_sum_products_gt(small, ns, large, nl, preL, sumL, (int)n, &cnt1, &sum1);
    count_sum_products_gt(small, ns, large, nl, preL, sumL, (int)(n / two), &cnt2, &sum2);

    i64 count_total = cnt1 + cnt2;
    i64 sum_d = sum1 + (i64)two * sum2;
    return sum_d - n * count_total;
}

i64 p834_native(void) {
    int N = 1234567;
    build_spf(N);

    i64 total = 0;
    int lastA = -1, lastB = -1;
    DivData dataA = {0}, dataB = {0};

    for (int n = 3; n <= N; n++) {
        int even, B;
        if (n & 1) {
            even = n - 1;
            B = n;
        } else {
            even = n;
            B = n - 1;
        }

        /* odd part and 2^v2 */
        int two = even & (-even);
        int A = even / two;

        if (A != lastA) {
            if (lastA != -1) free_divdata(&dataA);
            dataA = divisors_with_prefix(A);
            lastA = A;
        }
        if (B != lastB) {
            if (lastB != -1) free_divdata(&dataB);
            dataB = divisors_with_prefix(B);
            lastB = B;
        }

        total += compute_T_from_divdata((i64)n, two, &dataA, &dataB);
    }

    free_divdata(&dataA);
    free_divdata(&dataB);
    free(spf);
    return total;
}
