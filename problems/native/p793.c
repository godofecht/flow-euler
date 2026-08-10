#include <stdint.h>
#include <stdlib.h>
typedef long long i64;
static int cmp_i64(const void *a, const void *b) {
    i64 x = *(const i64 *)a, y = *(const i64 *)b;
    return (x > y) - (x < y);
}
static i64 count_at_most(i64 *v, int n, i64 thr) {
    i64 count = 0;
    int right = n - 1;
    for (int left = 0; left < n; left++) {
        while (right > left && v[left] * v[right] > thr) right--;
        if (right <= left) break;
        count += right - left;
    }
    return count;
}
long long pe_solve(void) {
    int n = 1000003;
    i64 *values = malloc((size_t)n * sizeof(i64));
    i64 s0 = 290797;
    for (int i = 0; i < n; i++) {
        values[i] = s0;
        s0 = s0 * s0 % 50515093;
    }
    qsort(values, (size_t)n, sizeof(i64), cmp_i64);
    i64 lo = values[0] * values[1];
    i64 hi = values[n - 1] * values[n - 2];
    i64 target = ((i64)n * (n - 1) / 2 + 1) / 2;
    while (lo < hi) {
        i64 mid = lo + (hi - lo) / 2;
        if (count_at_most(values, n, mid) >= target) hi = mid;
        else lo = mid + 1;
    }
    free(values);
    return lo;
}
