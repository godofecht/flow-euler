#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef int64_t i64;
enum { MOD = 123454321, TARGET_RANK = 1000000, TARGET_OMEGA = 1000000 };

static int *PRIMES; static int PC;
static double *PRIME_LOGS;
static double LOG3;

static i64 *values; static int vcount, vcap;
static i64 count_only;
static int collect_flag;
static i64 limit_i;
static double log_limit;
static int threshold_g;

static void search(int start_index, i64 product, int used, double log_product) {
    if (used >= threshold_g) {
        count_only++;
        if (collect_flag) {
            if (vcount >= vcap) {
                vcap = vcap ? vcap * 2 : 1024;
                values = realloc(values, vcap * sizeof(i64));
            }
            values[vcount++] = product;
        }
    }
    int remaining = threshold_g - used;
    for (int index = start_index; index < PC; index++) {
        i64 prime = PRIMES[index];
        if (product > limit_i / prime) break;
        i64 next_product = product * prime;
        if (remaining > 0 && log_product + remaining * PRIME_LOGS[index] > log_limit + 1e-12)
            break;
        search(index, next_product, used + 1, log_product + PRIME_LOGS[index]);
    }
}

static i64 values_count(int threshold, int collect) {
    threshold_g = threshold;
    limit_i = 1;
    for (int i = 0; i < threshold; i++) limit_i *= 3;
    log_limit = threshold * LOG3;
    collect_flag = collect;
    count_only = 0;
    vcount = 0;
    search(0, 1, 0, 0.0);
    return count_only;
}

static int cmp_i64(const void *a, const void *b) {
    i64 x = *(const i64*)a, y = *(const i64*)b;
    return (x > y) - (x < y);
}

static i64 mod_pow(i64 a, i64 e) {
    i64 r = 1; a %= MOD;
    while (e > 0) {
        if (e & 1) r = (i64)((__int128)r * a % MOD);
        a = (i64)((__int128)a * a % MOD);
        e >>= 1;
    }
    return r;
}

long long pe615_answer(void) {
    int lim = 300000;
    char *sieve = calloc(lim + 1, 1);
    PRIMES = calloc(lim, sizeof(int));
    PC = 0;
    for (int i = 2; i <= lim; i++) {
        if (!sieve[i]) {
            PRIMES[PC++] = i;
            if ((i64)i * i <= lim)
                for (int j = i * i; j <= lim; j += i) sieve[j] = 1;
        }
    }
    PRIME_LOGS = calloc(PC, sizeof(double));
    for (int i = 0; i < PC; i++) PRIME_LOGS[i] = log((double)PRIMES[i]);
    LOG3 = log(3.0);
    values = NULL; vcap = 0;

    int threshold = 1;
    i64 value = 0;
    while (1) {
        threshold++;
        i64 cnt = values_count(threshold, 0);
        if (cnt >= TARGET_RANK) {
            values_count(threshold, 1);
            qsort(values, vcount, sizeof(i64), cmp_i64);
            value = values[TARGET_RANK - 1];
            break;
        }
    }
    i64 ans = (i64)((__int128)(value % MOD) * mod_pow(2, TARGET_OMEGA - threshold) % MOD);
    free(sieve); free(PRIMES); free(PRIME_LOGS); free(values);
    return ans;
}
