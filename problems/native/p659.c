#include <stdint.h>
#include <stdlib.h>

typedef int64_t i64;
enum { LIMIT = 10000000 };

long long pe659_answer(void) {
    i64 *f = calloc(LIMIT + 1, sizeof(i64));
    i64 *maxelem = calloc(LIMIT + 1, sizeof(i64));
    if (!f || !maxelem) return 0;
    for (i64 x = 0; x <= LIMIT; x++) f[x] = 4 * x * x + 1;
    for (i64 x = 1; x <= LIMIT; x++) {
        i64 div = f[x];
        if (div > 1) {
            for (i64 curr = x % div; curr <= LIMIT; curr += div) {
                if (f[curr] % div == 0) {
                    if (maxelem[curr] < div) maxelem[curr] = div;
                    while (f[curr] % div == 0) f[curr] /= div;
                }
            }
            for (i64 curr = (-x % div + div) % div; curr <= LIMIT; curr += div) {
                if (curr == 0) continue;
                if (f[curr] % div == 0) {
                    if (maxelem[curr] < div) maxelem[curr] = div;
                    while (f[curr] % div == 0) f[curr] /= div;
                }
            }
        }
    }
    i64 sum = 0;
    const i64 MOD = 1000000000000000000LL;
    for (i64 i = 0; i <= LIMIT; i++) sum = (sum + maxelem[i]) % MOD;
    free(f); free(maxelem);
    return sum;
}
