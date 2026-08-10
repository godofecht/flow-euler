#include <stdint.h>
#include <stdlib.h>
#include <math.h>

typedef int64_t i64;
static i64 N;
static int *primes, pc;

static i64 generate(i64 curr, int idx) {
    i64 total = 0;
    i64 p = primes[idx];
    if (curr > N / (p * p * p)) return 0;
    curr *= p * p * p;
    while (curr <= N) {
        total += N / curr;
        for (int j = idx + 1; j < pc; j++) {
            i64 t = generate(curr, j);
            if (t == 0) break;
            total += t;
        }
        if (curr > N / p) break;
        curr *= p;
    }
    return total;
}

long long pe694_answer(void) {
    N = 1000000000000000000LL;
    i64 lim = (i64)(cbrt((double)N) + 2);
    char *sieve = calloc(lim + 1, 1);
    primes = calloc(lim, sizeof(int));
    pc = 0;
    for (i64 i = 2; i <= lim; i++) {
        if (!sieve[i]) {
            primes[pc++] = (int)i;
            for (i64 j = i * i; j <= lim; j += i) sieve[j] = 1;
        }
    }
    i64 ans = N;
    for (int i = 0; i < pc; i++) ans += generate(1, i);
    free(sieve); free(primes);
    return ans;
}
