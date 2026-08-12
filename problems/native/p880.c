// Project Euler 880: Nested Radicals
// Compute H(10^15) mod (1031^3 + 2).
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MOD (1031LL*1031LL*1031LL + 2)
#define N 1000000000000000LL // 10^15

static long long icbrt(long long n) {
    if (n <= 1) return n;
    long long r = (long long)round(pow((double)n, 1.0/3.0));
    while ((r+1)*(r+1)*(r+1) <= n) r++;
    while (r*r*r > n) r--;
    return r;
}

static long long iroot4(long long n) {
    long long r = (long long)sqrt((double)sqrt((double)n));
    while ((__int128)(r+1)*(r+1)*(r+1)*(r+1) <= n) r++;
    while ((__int128)r*r*r*r > n) r--;
    return r;
}

static long long isqrt_ll(long long n) {
    if (n < 0) return 0;
    long long r = (long long)sqrt((double)n);
    while ((r+1)*(r+1) <= n) r++;
    while (r*r > n) r--;
    return r;
}

static int gcd(int a, int b) {
    while (b) { int t = a % b; a = b; b = t; }
    return a;
}

static int *spf_array;

static void build_spf(int limit) {
    spf_array = malloc((limit + 1) * sizeof(int));
    for (int i = 0; i <= limit; i++) spf_array[i] = i;
    spf_array[1] = 1;
    for (int p = 2; (long long)p * p <= limit; p++) {
        if (spf_array[p] != p) continue;
        for (int m = p * p; m <= limit; m += p) {
            if (spf_array[m] == m) spf_array[m] = p;
        }
    }
}

static int *cube_free_table(int limit) {
    build_spf(limit);
    int *cf = malloc((limit + 1) * sizeof(int));
    cf[0] = 1; cf[1] = 1;
    for (int n = 2; n <= limit; n++) {
        int p = spf_array[n];
        int m = n / p;
        int e = 1;
        while (m % p == 0) { m /= p; e++; }
        int rem = e % 3;
        if (rem == 0) cf[n] = cf[m];
        else if (rem == 1) cf[n] = cf[m] * p;
        else cf[n] = cf[m] * p * p;
    }
    return cf;
}

static long long sumsq_mod(long long k, long long mod) {
    // k*(k+1)*(2k+1)/6 mod mod, using __int128 for exact computation
    __int128 val = (__int128)k * (k + 1) * (2 * k + 1) / 6;
    return (long long)(val % mod);
}

static long long H_mod(long long limit) {
    long long b_limit = iroot4(4 * limit);

    long long max_odd_a = 0;
    {
        long long c = icbrt(limit);
        if (c > 1) max_odd_a = (c - 1) / 4;
    }
    long long max_even_a = 0;
    {
        long long c = icbrt(limit / 4);
        if (c > 1) max_even_a = (c - 1) / 2;
    }
    long long cf_limit_val = 4 * max_odd_a;
    if (max_even_a > cf_limit_val) cf_limit_val = max_even_a;
    if (2 * b_limit > cf_limit_val) cf_limit_val = 2 * b_limit;

    int cf_limit = (int)cf_limit_val;
    int *cf = cube_free_table(cf_limit);

    int *cf4 = malloc((max_odd_a + 1) * sizeof(int));
    cf4[0] = 0;
    for (long long a = 1; a <= max_odd_a; a++)
        cf4[a] = cf[4 * a];

    long long total = 0;
    long long mod = MOD;

    // Odd b
    for (long long b = 1; b <= b_limit; b += 2) {
        long long a_limit = (icbrt(limit / b) - b) / 4;
        int cf_b = cf[b];
        for (long long a = 1; a <= a_limit; a++) {
            if (gcd((int)a, (int)b) != 1 || cf_b == cf4[a])
                continue;

            long long x_base = b + 4 * a;
            __int128 x = (__int128)b * x_base * x_base * x_base;
            long long y_base = a - 2 * b;
            __int128 yc = (__int128)y_base * y_base * y_base;
            __int128 y_abs = 4 * (__int128)a * (yc < 0 ? -yc : yc);

            __int128 max_coord = (x >= y_abs) ? x : y_abs;
            if (max_coord > (__int128)limit) continue;

            long long mc = (long long)max_coord;
            long long tmax = isqrt_ll(limit / mc);

            long long s = sumsq_mod(tmax, mod);
            long long xy_mod = (long long)((x + y_abs) % (__int128)mod);
            long long contribution = (xy_mod * s) % mod;
            total = (total + contribution) % mod;
        }
    }

    // Even b
    for (long long b = 2; b <= b_limit; b += 2) {
        long long half_b = b / 2;
        long long a_limit = (icbrt(limit / (2 * b)) - half_b) / 2;
        int cf_2b = cf[2 * b];
        for (long long a = 1; a <= a_limit; a += 2) {
            if (gcd((int)a, (int)b) != 1 || cf_2b == cf[a])
                continue;

            long long x_base = half_b + 2 * a;
            __int128 x = (__int128)2 * b * x_base * x_base * x_base;
            long long y_base = a - 2 * b;
            __int128 yc = (__int128)y_base * y_base * y_base;
            __int128 y_abs = (__int128)a * (yc < 0 ? -yc : yc);

            __int128 max_coord = (x >= y_abs) ? x : y_abs;
            if (max_coord > (__int128)limit) continue;

            long long mc = (long long)max_coord;
            long long tmax = isqrt_ll(limit / mc);

            long long s = sumsq_mod(tmax, mod);
            long long xy_mod = (long long)((x + y_abs) % (__int128)mod);
            long long contribution = (xy_mod * s) % mod;
            total = (total + contribution) % mod;
        }
    }

    free(cf);
    free(cf4);
    free(spf_array);
    return total;
}

long long p880_native(void) {
    return H_mod(N) % MOD;
}
