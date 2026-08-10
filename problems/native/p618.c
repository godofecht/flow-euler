#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;
enum { LIMIT = 46368, MOD = 1000000000LL };

long long pe618_answer(void) {
    char *prime = calloc(LIMIT + 1, 1);
    i64 *d = calloc(LIMIT + 1, sizeof(i64));
    if (!prime || !d) return 0;
    for (int i = 0; i <= LIMIT; i++) prime[i] = 1;
    prime[0] = prime[1] = 0;
    for (int i = 2; i * (i64)i <= LIMIT; i++) if (prime[i]) {
        for (int j = i * i; j <= LIMIT; j += i) prime[j] = 0;
    }
    d[0] = 1;
    for (int p = 2; p <= LIMIT; p++) if (prime[p]) {
        for (int i = p; i <= LIMIT; i++) {
            d[i] = (d[i] + (i64)p * d[i - p]) % MOD;
        }
    }
    i64 total = 0;
    i64 a = 1, b = 2;
    while (b <= LIMIT) {
        total = (total + d[b]) % MOD;
        i64 c = a + b; a = b; b = c;
    }
    free(prime); free(d);
    return total;
}
