/* Fair C twin: same algorithm + in-process multi-round timing. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

static int64_t sum_primes_below(int64_t limit) {
    char *sieve = calloc((size_t)limit, 1);
    if (!sieve) return -1;
    sieve[0] = 1;
    sieve[1] = 1;
    for (int64_t p = 2; p * p < limit; p++) {
        if (sieve[p] == 0) {
            for (int64_t m = p * p; m < limit; m += p)
                sieve[m] = 1;
        }
    }
    int64_t total = 0;
    for (int64_t i = 2; i < limit; i++)
        if (sieve[i] == 0) total += i;
    free(sieve);
    return total;
}

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

int main(int argc, char **argv) {
    int rounds = 200;
    if (argc > 1) rounds = atoi(argv[1]);
    const int64_t limit = 2000000;
    volatile int64_t sink = sum_primes_below(limit);
    double t0 = now_s();
    for (int i = 0; i < rounds; i++)
        sink = sum_primes_below(limit);
    double t1 = now_s();
    printf("%lld\n", (long long)sink);
    fprintf(stderr, "c: %.4fs / %d rounds (%.3f ms/run)\n",
            t1 - t0, rounds, 1000.0 * (t1 - t0) / rounds);
    return 0;
}
