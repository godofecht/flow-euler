#include <stdint.h>
#include <stdlib.h>
#include <math.h>

typedef long long i64;

static int cmp_i64(const void *a, const void *b) {
    i64 x = *(const i64 *)a, y = *(const i64 *)b;
    return (x > y) - (x < y);
}

long long p757_native(void) {
    i64 limit = 100000000000000LL; /* 10^14 */
    size_t cap = 1 << 24; /* 16M to start, grows as needed */
    i64 *arr = (i64 *)malloc(cap * sizeof(i64));
    size_t cnt = 0;

    for (i64 x = 1;; x++) {
        i64 xx = x * (x + 1);
        if (xx * xx > limit) break; /* x <= y, so smallest product is xx*xx */
        i64 rem = limit / xx;
        /* y(y+1) <= rem  =>  y <= (sqrt(4*rem+1)-1)/2 */
        i64 ymax = (i64)((sqrt((double)(4 * rem + 1)) - 1.0) / 2.0);
        while (ymax * (ymax + 1) > rem) ymax--;
        for (i64 y = x; y <= ymax; y++) {
            i64 n = xx * y * (y + 1);
            if (cnt >= cap) {
                cap <<= 1;
                arr = (i64 *)realloc(arr, cap * sizeof(i64));
            }
            arr[cnt++] = n;
        }
    }

    qsort(arr, cnt, sizeof(i64), cmp_i64);

    i64 distinct = 0;
    for (size_t i = 0; i < cnt;) {
        distinct++;
        size_t j = i + 1;
        while (j < cnt && arr[j] == arr[i]) j++;
        i = j;
    }

    free(arr);
    return distinct;
}
