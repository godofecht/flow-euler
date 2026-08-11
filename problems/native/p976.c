// Project Euler 976: F(10^7, 10^7) modulo 1234567891.
// Port of the reference Python solver.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MOD 1234567891LL

static long long *build_inverses(long long n, long long mod) {
    long long *inv = (long long *)malloc((size_t)(n + 1) * sizeof(long long));
    if (!inv) {
        fprintf(stderr, "oom\n");
        exit(1);
    }
    memset(inv, 0, (size_t)(n + 1) * sizeof(long long));
    if (n >= 1) inv[1] = 1;
    for (long long i = 2; i <= n; i++) {
        inv[i] = mod - (mod / i) * inv[mod % i] % mod;
    }
    return inv;
}

long long p976_native(void) {
    long long n = 10000000LL;
    long long k = 10000000LL;

    long long e = n / 2;           // 5000000
    long long a_cnt = (n + 3) / 4; // 2500000
    long long b_cnt = (n + 1) / 4; // 2500000
    long long c = b_cnt - a_cnt;   // 0

    if (e == 0) {
        // Only odd strips exist. Not reached for n = 10^7, but kept for parity
        // with the reference solver.
        long long *inv = build_inverses(k + 2, MOD);
        long long inv2 = (MOD + 1) / 2;
        long long h = 1, qv = 1, sum_odd_a = 0, ans = 0;
        for (long long s = 0; s <= k; s++) {
            if (s > 0) {
                h = h * (a_cnt + b_cnt + s - 1) % MOD * inv[s] % MOD;
                if (s % 2 == 0) {
                    long long r = s / 2;
                    qv = qv * (a_cnt + r - 1) % MOD * inv[r] % MOD;
                }
            }
            long long coeff;
            if (c == 0) {
                coeff = (s % 2 == 0) ? qv : 0;
            } else if (c == 1) {
                coeff = qv;
            } else {
                coeff = (s % 2 == 0) ? qv : (MOD - qv);
            }
            long long diff = (h - coeff) % MOD;
            if (diff < 0) diff += MOD;
            long long h_odd_a = diff * inv2 % MOD;
            sum_odd_a = (sum_odd_a + h_odd_a) % MOD;
            ans = sum_odd_a;
        }
        free(inv);
        return ans % MOD;
    }

    long long max_inv = e + k + 2;
    long long *inv = build_inverses(max_inv, MOD);
    long long inv2 = (MOD + 1) / 2;

    // total_even at m = k: C(e + k - 1, k)
    long long total_even = 1;
    for (long long m = 0; m < k; m++) {
        total_even = total_even * (e + m) % MOD * inv[m + 1] % MOD;
    }

    // e0 at q = k / 2: C(e + q - 1, q)
    long long qmax = k / 2;
    long long e0 = 1;
    for (long long q = 0; q < qmax; q++) {
        e0 = e0 * (e + q) % MOD * inv[q + 1] % MOD;
    }

    long long h = 1;     // h_s
    long long qv = 1;    // C(a_cnt + r - 1, r)
    long long qsum = 1;  // C(a_cnt + r, r)
    long long sum_even = 0, sum_odd = 0, sum_odd_a = 0, ans = 0;
    long long ab = a_cnt + b_cnt;

    for (long long s = 0; s <= k; s++) {
        if (s > 0) {
            h = h * (ab + s - 1) % MOD * inv[s] % MOD;
            if (s % 2 == 0) {
                long long r = s / 2;
                qv = qv * (a_cnt + r - 1) % MOD * inv[r] % MOD;
                if (c == 1) {
                    qsum = qsum * (a_cnt + r) % MOD * inv[r] % MOD;
                }
            }
        }

        long long coeff;
        if (c == 0) {
            coeff = (s % 2 == 0) ? qv : 0;
        } else if (c == 1) {
            coeff = qsum;
        } else {
            coeff = (s % 2 == 0) ? qv : (MOD - qv);
        }

        long long diff = (h - coeff) % MOD;
        if (diff < 0) diff += MOD;
        long long h_odd_a = diff * inv2 % MOD;

        if (s % 2 == 0) {
            sum_even = (sum_even + h) % MOD;
        } else {
            sum_odd = (sum_odd + h) % MOD;
        }
        sum_odd_a = (sum_odd_a + h_odd_a) % MOD;

        long long m = k - s;
        long long t0, t1;
        if (m % 2 == 0) {
            long long e0_m = e0;
            t0 = e0_m * sum_odd_a % MOD;
            long long d1 = (total_even - e0_m) % MOD;
            if (d1 < 0) d1 += MOD;
            t1 = d1 * sum_even % MOD;
        } else {
            t0 = 0;
            t1 = total_even * sum_odd % MOD;
        }
        ans = (ans + t0 + t1) % MOD;

        if (m > 0) {
            total_even = total_even * m % MOD * inv[e + m - 1] % MOD;
        }
        if (m % 2 == 0 && m >= 2) {
            long long qcur = m / 2;
            e0 = e0 * qcur % MOD * inv[e + qcur - 1] % MOD;
        }
    }

    free(inv);
    return ans % MOD;
}
