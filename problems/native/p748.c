/*
 * Project Euler 748 - Upside Down Diophantine Equation
 *
 * Compute the last 9 digits of S(10^16) where S(N) = sum(x+y+z)
 * over primitive integer solutions of 1/x^2 + 1/y^2 = 13/z^2
 * with 1 <= x,y,z <= N and x <= y.
 *
 * Algorithm: parametrize via Gaussian integers.
 *   (p + iq) = (3 + 2i) * (m + in)^2,  r = m^2 + n^2,
 *   then x = q*r, y = p*r, z = p*q.
 * Enumerate coprime (m,n) with opposite parity, skip the
 * non-primitive family where 13 | p and 13 | q.
 */

#include <stdint.h>
#include <stdlib.h>

typedef unsigned __int128 u128;

/* GCD of non-negative inputs. */
static long long gcd_ll(long long a, long long b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) {
        long long t = a % b;
        a = b;
        b = t;
    }
    return a;
}

/* Floor of integer square root for 128-bit values. */
static u128 isqrt128(u128 n) {
    if (n == 0) return 0;
    u128 x = n;
    u128 y = (x + 1) / 2;
    while (y < x) {
        x = y;
        y = (x + n / x) / 2;
    }
    return x;
}

/* Floor of fourth root of a 128-bit value. */
static u128 fourth_root_floor128(u128 n) {
    if (n == 0) return 0;
    u128 x = isqrt128(isqrt128(n));
    u128 xp1 = x + 1;
    while (xp1 * xp1 * xp1 * xp1 <= n) {
        x = xp1;
        xp1 = x + 1;
    }
    while (x * x * x * x > n) x--;
    return x;
}

/* Floor of integer square root for 64-bit values. */
static long long isqrt64(long long n) {
    if (n <= 0) return 0;
    long long x = n;
    long long y = (x + 1) / 2;
    while (y < x) {
        x = y;
        y = (x + n / x) / 2;
    }
    return x;
}

long long p748_native(void) {
    long long N = 10000000000000000LL; /* 10^16 */
    long long mod = 1000000000LL;      /* 10^9  */

    /* r_max = floor((2*N^2 / 13)^(1/4)) */
    u128 Nsq = (u128)N * (u128)N;
    u128 val = (u128)2 * Nsq / (u128)13;
    long long r_max = (long long)fourth_root_floor128(val);

    long long m_max = isqrt64(r_max);

    long long total = 0;
    /* Reduce before i64 overflow: each s <= ~3*10^16. */
    long long THRESH = 1000000000000000000LL; /* 10^18 */

    for (long long m = 1; m <= m_max; m++) {
        long long mm = m * m;
        long long n_max = isqrt64(r_max - mm);

        /* Opposite parity: m and n must differ in parity. */
        long long n_start = (m & 1) ? 0 : 1;

        for (long long n = n_start; n <= n_max; n += 2) {
            if (gcd_ll(m, n) != 1) continue;

            long long nn = n * n;
            long long r = mm + nn;

            long long u = mm - nn;   /* m^2 - n^2 */
            long long v = 2 * m * n; /* 2mn        */

            /* (3 + 2i)(u + iv) = (3u - 2v) + i(3v + 2u) */
            long long a = 3 * u - 2 * v;
            if (a < 0) a = -a;
            long long b = 3 * v + 2 * u;
            if (b < 0) b = -b;

            /* Order so p >= q (gives x <= y after mapping). */
            long long p, q;
            if (a < b) { p = b; q = a; }
            else       { p = a; q = b; }

            /* Non-primitive family: 13 | p and 13 | q. */
            if (p % 13 == 0 && q % 13 == 0) continue;

            long long x = q * r;
            long long y = p * r;
            if (x > N || y > N) continue;
            long long z = p * q;
            if (z > N) continue;

            total += x + y + z;
            if (total >= THRESH) total %= mod;
        }
    }

    return total % mod;
}
