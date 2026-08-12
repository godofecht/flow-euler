/* Project Euler 809: Rational Recurrence Relation
   f(22/7) mod 10^15 via Ackermann-Peter function and tetration.
   Port of the Python reference solver. */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef __int128 i128;

/* Extended GCD: returns g, sets x and y such that a*x + b*y = g */
static long long egcd(long long a, long long b, long long *x, long long *y) {
    if (b == 0) { *x = 1; *y = 0; return a; }
    long long x1, y1;
    long long g = egcd(b, a % b, &x1, &y1);
    *x = y1;
    *y = x1 - (a / b) * y1;
    return g;
}

static long long inv_mod(long long a, long long m) {
    long long x, y;
    egcd(a % m, m, &x, &y);
    if (x < 0) x += m;
    return x % m;
}

static long long crt(long long r1, long long m1, long long r2, long long m2) {
    /* x ≡ r1 (mod m1), x ≡ r2 (mod m2), m1 and m2 coprime */
    i128 diff = (i128)((r2 - r1) % m2);
    if (diff < 0) diff += m2;
    i128 k = diff * inv_mod(m1 % m2, m2) % m2;
    i128 result = (i128)r1 + (i128)m1 * k;
    /* result < m1 * m2 */
    return (long long)(result % ((i128)m1 * m2));
}

/* factor n = 2^a * 5^b * rem */
static void factor_2_5(long long n, int *a, int *b, long long *rem) {
    *a = 0;
    while (n % 2 == 0) { (*a)++; n /= 2; }
    *b = 0;
    while (n % 5 == 0) { (*b)++; n /= 5; }
    *rem = n;
}

static long long phi_2_5(long long n) {
    int a, b;
    long long rem;
    factor_2_5(n, &a, &b, &rem);
    long long res = n;
    if (a) res /= 2;
    if (b) res = (res / 5) * 4;
    return res;
}

static int totient_chain_len(long long n) {
    int steps = 0;
    while (n != 1) {
        n = phi_2_5(n);
        steps++;
    }
    return steps;
}

/* min(2↑↑height, cap) for base=2 */
static long long tetration_cap(int height, long long cap) {
    if (cap <= 0) return 0;
    long long v = 2;
    if (height <= 1) return (v < cap) ? v : cap;
    for (int i = 2; i <= height; i++) {
        if (v >= 60) return cap;
        v = 1LL << v;
        if (v >= cap) return cap;
    }
    return v;
}

/* 2↑↑height (mod 2^a) */
static long long tetration_mod_pow2(int height, int a) {
    if (a <= 0) return 0;
    long long mod = 1LL << a;
    if (height == 1) return 2 % mod;
    long long exp = tetration_cap(height - 1, a);
    if (exp >= a) return 0;
    return (1LL << exp) % mod;
}

static long long powmod_128(long long base, long long exp, long long mod) {
    i128 r = 1 % mod;
    i128 b = base % mod;
    while (exp) {
        if (exp & 1) r = r * b % mod;
        b = b * b % mod;
        exp >>= 1;
    }
    return (long long)r;
}

/* 2↑↑height (mod mod) for mod of the form 2^a * 5^b */
static long long tetration_mod(int height, long long mod) {
    if (mod == 1) return 0;
    int a, b;
    long long rem;
    factor_2_5(mod, &a, &b, &rem);

    if (b == 0) return tetration_mod_pow2(height, a);

    if (a == 0) {
        /* Odd modulus (5^b) */
        if (height == 1) return 2 % mod;
        long long exp = tetration_mod(height - 1, phi_2_5(mod));
        return powmod_128(2, exp, mod);
    }

    /* Mixed modulus: CRT */
    long long m2 = 1LL << a;
    long long m5 = 1;
    for (int i = 0; i < b; i++) m5 *= 5;
    long long r2 = tetration_mod_pow2(height, a);
    long long r5 = tetration_mod(height, m5);
    return crt(r2, m2, r5, m5) % (m2 * m5);
}

static long long stable_tetration_mod(long long mod) {
    int height = totient_chain_len(mod) + 1;
    return tetration_mod(height, mod);
}

long long p809_native(void) {
    long long mod2 = 1LL << 15;
    long long mod5 = 1;
    for (int i = 0; i < 15; i++) mod5 *= 5; /* 5^15 = 30517578125 */

    /* 2 ↑↑↑↑ 6 is a power of two, so mod 2^15 it is 0, hence A(6,3) ≡ -3 */
    long long r2 = (-3) % mod2;
    if (r2 < 0) r2 += mod2;

    /* For the odd part: tall tetration of 2's modulo 5^15 */
    long long tower_mod5 = stable_tetration_mod(mod5);
    long long r5 = (tower_mod5 - 3) % mod5;
    if (r5 < 0) r5 += mod5;

    /* CRT to combine mod 2^15 and mod 5^15 → mod 10^15 */
    long long result = crt(r2, mod2, r5, mod5);

    /* result mod 10^15 */
    long long mod10_15 = mod2 * mod5; /* 10^15 */
    return result % mod10_15;
}
