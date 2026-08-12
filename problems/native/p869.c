/* Project Euler 869: Prime Guessing
 *
 * E(10^8): expected answer 14.97696693
 *
 * Build a binary trie keyed by the reversed binary representation of each
 * prime (LSB first).  At every internal node the score contribution is
 * max(count_left, count_right).  The answer is total_score / pi(N).
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef long long i64;

static i64 total_points = 0;
static int *tmpbuf;

static void calc_points(int *arr, int lo, int hi, int bit) {
    if (hi <= lo) return;

    int count0 = 0, count1 = 0;
    for (int i = lo; i < hi; i++) {
        int v = arr[i] >> bit;
        if (v == 0) continue;
        if (v & 1) count1++;
        else count0++;
    }

    if (count0 + count1 == 0) return;

    total_points += (i64)((count0 > count1) ? count0 : count1);

    int p0 = 0, p1 = count0, p2 = count0 + count1;
    for (int i = lo; i < hi; i++) {
        int v = arr[i] >> bit;
        if (v == 0) tmpbuf[p2++] = arr[i];
        else if (v & 1) tmpbuf[p1++] = arr[i];
        else tmpbuf[p0++] = arr[i];
    }

    memcpy(arr + lo, tmpbuf, (size_t)(hi - lo) * sizeof(int));

    int mid0 = lo + count0;
    int mid1 = mid0 + count1;
    calc_points(arr, lo, mid0, bit + 1);
    calc_points(arr, mid0, mid1, bit + 1);
}

double p869_native(void) {
    int N = 100000000;

    char *is_prime = malloc((size_t)(N + 1));
    memset(is_prime, 1, (size_t)(N + 1));
    is_prime[0] = is_prime[1] = 0;
    int limit = (int)sqrt((double)N);
    for (int i = 2; i <= limit; i++) {
        if (is_prime[i]) {
            for (int j = i * i; j <= N; j += i) is_prime[j] = 0;
        }
    }

    int pc = 0;
    for (int i = 2; i <= N; i++) if (is_prime[i]) pc++;

    int *primes = malloc((size_t)pc * sizeof(int));
    int idx = 0;
    for (int i = 2; i <= N; i++) if (is_prime[i]) primes[idx++] = i;
    free(is_prime);

    tmpbuf = malloc((size_t)pc * sizeof(int));

    total_points = 0;
    calc_points(primes, 0, pc, 0);

    double result = (double)total_points / (double)pc;

    free(primes);
    free(tmpbuf);

    return result;
}
