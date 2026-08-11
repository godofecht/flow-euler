#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef long long i64;

/* Count k in [0..n^2] achievable with complexity <= 3. */
static i64 count_le_3(int n) {
    i64 n2 = (i64)n * n;
    i64 half = n2 / 2;
    unsigned char *seen = calloc((size_t)(half + 1), 1);

    /* Orbit 1: k = x*y, 0 <= x <= y <= n. */
    for (int y = 1; y <= n; y++) {
        for (int x = 0; x <= y; x++) {
            i64 k = (i64)x * y;
            if (k > half) k = n2 - k;
            seen[k] = 1;
        }
    }
    seen[0] = 1;

    /* Orbit 2: k = v^2 - d^2, 0 <= d <= v <= n. */
    for (int v = 0; v <= n; v++) {
        i64 vv = (i64)v * v;
        for (int d = 0; d <= v; d++) {
            i64 k = vv - (i64)d * d;
            if (k > half) k = n2 - k;
            if (k >= 0 && k <= half) seen[k] = 1;
        }
    }

    /* q[b] = 2*b*(n-b) for b = 0..floor(n/2). */
    int bmax_global = n / 2;
    i64 *q = malloc((size_t)(bmax_global + 1) * sizeof(i64));
    {
        i64 cur = 0, delta = 2 * (i64)(n - 1);
        for (int b = 0; b < bmax_global; b++) {
            q[b] = cur;
            cur += delta;
            delta -= 4;
        }
        q[bmax_global] = cur;
    }

    /* Orbits 3/4/5/6: depend on s = a+b and ab = a(s-a). */
    for (int s = 0; s <= n; s++) {
        int c = n - s;
        i64 c2 = (i64)c * c;
        i64 c2n = (i64)c * (2 * n - c);
        i64 cs = (i64)c * s;

        if ((s & 1) == 0) {
            int v = s / 2;
            i64 base = (i64)v * v;
            for (int d = 0; d <= v; d++) {
                i64 ab = base - (i64)d * d;
                i64 two_ab = ab << 1;
                i64 k;
                k = c2 + two_ab;  if (k > half) k = n2 - k; seen[k] = 1;
                k = c2n + two_ab; if (k > half) k = n2 - k; seen[k] = 1;
                k = cs + ab;      if (k > half) k = n2 - k; seen[k] = 1;
                k = (cs + ab) << 1; if (k > half) k = n2 - k; seen[k] = 1;
            }
        } else {
            int v = s / 2;
            i64 base = (i64)v * (v + 1);
            i64 pr = 0, step = 2;
            for (int d = 0; d <= v; d++) {
                i64 ab = base - pr;
                i64 two_ab = ab << 1;
                i64 k;
                k = c2 + two_ab;  if (k > half) k = n2 - k; seen[k] = 1;
                k = c2n + two_ab; if (k > half) k = n2 - k; seen[k] = 1;
                k = cs + ab;      if (k > half) k = n2 - k; seen[k] = 1;
                k = (cs + ab) << 1; if (k > half) k = n2 - k; seen[k] = 1;
                pr += step;
                step += 2;
            }
        }
    }

    /* Orbit 7: k = c^2 + 2*b*(n-b), 0 <= b <= min(n-c, n/2). */
    for (int c = 0; c <= n; c++) {
        i64 c2 = (i64)c * c;
        int bmax = n - c;
        if (bmax > bmax_global) bmax = bmax_global;
        for (int b = 0; b <= bmax; b++) {
            i64 k = c2 + q[b];
            if (k > half) k = n2 - k;
            seen[k] = 1;
        }
    }

    i64 cnt = 0;
    for (i64 i = 0; i <= half; i++) cnt += seen[i];
    i64 result = (n2 & 1) ? 2 * cnt : 2 * cnt - seen[half];
    free(seen);
    free(q);
    return result;
}

static int cmp_i64(const void *a, const void *b) {
    i64 x = *(const i64 *)a, y = *(const i64 *)b;
    return (x > y) - (x < y);
}

/* Count k in [0..n^2] whose minimum complexity is exactly 2. */
static i64 count_eq_2(int n) {
    i64 n2 = (i64)n * n;
    i64 *vals = malloc((size_t)(4 * (n + 1)) * sizeof(i64));
    int m = 0;
    for (int a = 0; a <= n; a++) {
        i64 a2 = (i64)a * a;
        i64 t = 2 * (i64)a * (n - a);
        vals[m++] = a2;
        vals[m++] = n2 - a2;
        vals[m++] = t;
        vals[m++] = n2 - t;
    }
    qsort(vals, m, sizeof(i64), cmp_i64);
    i64 cnt = 0;
    for (int i = 0; i < m; i++) {
        if (vals[i] == 0 || vals[i] == n2) continue;
        if (i == 0 || vals[i] != vals[i - 1]) cnt++;
    }
    free(vals);
    return cnt;
}

long long p782_native(void) {
    int n = 10000;
    i64 n2 = (i64)n * n;
    i64 total = n2 + 1;
    i64 s3 = count_le_3(n);
    i64 n4 = total - s3;
    i64 n2cnt = count_eq_2(n);
    return 3 * total - 4 - n2cnt + n4;
}
