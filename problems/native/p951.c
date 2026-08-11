#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef long long i64;

static i64 comb(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    i64 result = 1;
    for (int i = 0; i < k; i++) {
        result = result * (i64)(n - i) / (i64)(i + 1);
    }
    return result;
}

/* Count length 2n strings with n R, n B and no maximal run of length 2.
   dp[r][last][rs]: r = R used (0..n), last = 0/1, rs = run-length category
   1, 2, or 3 (>=3). Switching colours is forbidden when rs == 2.
   Final state must also avoid rs == 2. */
static i64 count_balanced_no_run2(int n) {
    if (n <= 0) return (n == 0) ? 1 : 0;
    int L = 2 * n;
    int R = n + 1;
    size_t sz = (size_t)R * 2 * 3;
    i64 *dp = calloc(sz, sizeof(i64));
    i64 *nw = calloc(sz, sizeof(i64));
#define IDX(r, last, rs) ((size_t)(r) * 2 * 3 + (last) * 3 + (rs))
    /* first card */
    dp[IDX(1, 0, 0)] = 1; /* start with R, run category 1 */
    dp[IDX(0, 1, 0)] = 1; /* start with B, run category 1 */

    for (int pos = 1; pos < L; pos++) {
        memset(nw, 0, sz * sizeof(i64));
        for (int r = 0; r <= n; r++) {
            int b = pos - r;
            if (b < 0 || b > n) continue;
            for (int last = 0; last < 2; last++) {
                for (int rs = 0; rs < 3; rs++) {
                    i64 val = dp[IDX(r, last, rs)];
                    if (!val) continue;
                    int rscat = rs + 1; /* 1, 2, 3 */

                    /* add an R */
                    if (r + 1 <= n) {
                        if (last == 0) {
                            int nrs = (rscat < 3) ? rscat + 1 : 3;
                            nw[IDX(r + 1, 0, nrs - 1)] += val;
                        } else if (rscat != 2) {
                            nw[IDX(r + 1, 0, 0)] += val;
                        }
                    }
                    /* add a B */
                    if (b + 1 <= n) {
                        if (last == 1) {
                            int nrs = (rscat < 3) ? rscat + 1 : 3;
                            nw[IDX(r, 1, nrs - 1)] += val;
                        } else if (rscat != 2) {
                            nw[IDX(r, 1, 0)] += val;
                        }
                    }
                }
            }
        }
        i64 *tmp = dp; dp = nw; nw = tmp;
    }

    /* accept balanced states (r == n) where final run is not length 2 */
    i64 total = 0;
    for (int last = 0; last < 2; last++) {
        total += dp[IDX(n, last, 0)]; /* rs == 1 */
        total += dp[IDX(n, last, 2)]; /* rs == 3 */
    }
    free(dp);
    free(nw);
    return total;
#undef IDX
}

long long p951_native(void) {
    int n = 26;
    i64 total = comb(2 * n, n);
    i64 unfair = count_balanced_no_run2(n);
    return total - unfair;
}
