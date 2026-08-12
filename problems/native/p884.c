// Project Euler 884: Removing Cubes
// Compute S(10^17) where D(n) is the number of steps obtained by repeatedly
// subtracting the largest perfect cube not exceeding the current value.
// Uses a cube-interval recurrence with memoization.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static long long icbrt(__int128 n) {
    if (n <= 0) return 0;
    double d = (double)n;
    long long x = (long long)round(pow(d, 1.0/3.0));
    // Adjust
    while ((__int128)(x+1)*(x+1)*(x+1) <= n) x++;
    while ((__int128)x*x*x > n) x--;
    return x;
}

static long long d_greedy(long long n) {
    long long steps = 0;
    while (n > 0) {
        long long k = icbrt(n);
        n -= k * k * k;
        steps++;
    }
    return steps;
}

// Base prefix: F(n) = sum_{m=0}^{n-1} D(m) for n in [0..limit]
static long long *base_pref;
static int base_pref_len;

static void build_base_prefix(int limit) {
    long long *d = calloc(limit, sizeof(long long));
    for (int i = 0; i < limit; i++)
        d[i] = d_greedy(i);
    base_pref = calloc(limit + 1, sizeof(long long));
    base_pref_len = limit + 1;
    long long s = 0;
    for (int i = 0; i < limit; i++) {
        s += d[i];
        base_pref[i + 1] = s;
    }
    free(d);
}

// Memoization for F
#define MEMO_LIMIT 5000000
static long long *memo;
static char *memo_valid;

static long long *delta_prefix;
static int delta_prefix_len;

static long long F(long long N) {
    if (N <= 0) return 0;
    if (N < base_pref_len)
        return base_pref[N];
    if (N <= MEMO_LIMIT) {
        int idx = (int)N;
        if (memo_valid[idx])
            return memo[idx];
    }

    long long K = icbrt(N - 1);
    __int128 K3 = (__int128)K * K * K;
    long long L = (long long)(N - K3);

    long long res = (N - 1) + delta_prefix[K - 1] + F(L);

    if (N <= MEMO_LIMIT) {
        int idx = (int)N;
        memo[idx] = res;
        memo_valid[idx] = 1;
    }
    return res;
}

static long long compute_S(long long N) {
    if (N <= 1) return 0;

    long long K_target = icbrt(N - 1);

    // Largest delta among k <= K_target
    __int128 max_delta = 3 * (__int128)K_target * K_target + 3 * (__int128)K_target + 1;
    long long K0_max = icbrt(max_delta - 1);

    delta_prefix = calloc(K0_max + 1, sizeof(long long));
    delta_prefix_len = K0_max + 1;

    build_base_prefix(64);

    memo = calloc(MEMO_LIMIT + 1, sizeof(long long));
    memo_valid = calloc(MEMO_LIMIT + 1, 1);

    long long total_F_delta = 0;
    for (long long k = 1; k < K_target; k++) {
        long long delta_k = 3 * k * k + 3 * k + 1;
        long long f_delta = F(delta_k);
        total_F_delta += f_delta;
        if (k <= K0_max)
            delta_prefix[k] = total_F_delta;
    }

    __int128 K3 = (__int128)K_target * K_target * K_target;
    long long L = (long long)(N - K3);
    long long result = (N - 1) + total_F_delta + F(L);

    free(base_pref);
    free(delta_prefix);
    free(memo);
    free(memo_valid);
    return result;
}

long long p884_native(void) {
    return compute_S(100000000000000000LL); // 10^17
}
