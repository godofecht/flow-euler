#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Project Euler 830 - Binomials and Powers.
 *
 * S(n) = sum_{k=0..n} C(n,k) * k^n, computed mod 83^3 * 89^3 * 97^3.
 * For n = 10^18, we use the identity:
 *   S(n) = sum_{j=0..J} C(n,j) * (j! * S(n,j)) * 2^(n-j)  mod p^a
 * where J is bounded by the p-adic valuation of the falling factorial.
 *
 * CRT combines the three prime-power results.
 */

typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

static i64 mmul(i64 a, i64 b, i64 mod) {
    return (i64)((i128)a * b % mod);
}

static i64 mpow(i64 base, i64 exp, i64 mod) {
    i64 result = 1 % mod;
    base %= mod;
    if (base < 0) base += mod;
    while (exp > 0) {
        if (exp & 1) result = mmul(result, base, mod);
        base = mmul(base, base, mod);
        exp >>= 1;
    }
    return result;
}

static i64 egcd(i64 a, i64 b, i64 *x, i64 *y) {
    if (b == 0) { *x = 1; *y = 0; return a; }
    i64 x1, y1;
    i64 g = egcd(b, a % b, &x1, &y1);
    *x = y1;
    *y = x1 - (a / b) * y1;
    return g;
}

static i64 modinv(i64 a, i64 mod) {
    i64 x, y;
    a %= mod;
    if (a < 0) a += mod;
    egcd(a, mod, &x, &y);
    return ((x % mod) + mod) % mod;
}

static int v_p(i64 x, int p) {
    if (x == 0) return 100;
    int c = 0;
    while (x % p == 0) { x /= p; c++; }
    return c;
}

/* Largest J such that v_p(n*(n-1)*...*(n-J+1)) < a */
static int max_j(i64 n, int p, int a) {
    if (n <= 0) return 0;
    int vp = 0;
    int j = 0;
    while (1) {
        if (vp >= a) return j - 1;
        j++;
        if ((i64)j > n) return j - 1;
        vp += v_p(n - j + 1, p);
    }
}

/* Forward differences of powers: F[j] = Delta^j f(0) for f(i) = i^n, mod p^a */
static void forward_diff(i64 n, int J, i64 mod, int p, i64 *F) {
    i64 *arr = (i64 *)malloc((size_t)(J + 1) * sizeof(i64));
    /* arr[i] = i^n mod p^a. For i divisible by p and n >= a, this is 0. */
    arr[0] = 0; /* n != 0, so 0^n = 0 */
    for (int i = 1; i <= J; i++) {
        if (i % p == 0) {
            arr[i] = 0; /* v_p(i^n) = n * v_p(i) >= n >> a */
        } else {
            arr[i] = mpow(i, n, mod);
        }
    }
    for (int j = 0; j <= J; j++) {
        F[j] = arr[0] % mod;
        for (int i = 0; i < J - j; i++) {
            arr[i] = (arr[i + 1] - arr[i]) % mod;
            if (arr[i] < 0) arr[i] += mod;
        }
    }
    free(arr);
}

/* C(n, 0..J) mod p^a using p-adic valuation tracking */
static void binom_prefix(i64 n, int J, i64 mod, int p, int a, i64 *out) {
    out[0] = 1 % mod;
    i64 u_c = 1 % mod;  /* p-free part of C(n,j) mod p^a */
    int v_c = 0;        /* p-adic valuation of C(n,j) */

    for (int j = 1; j <= J; j++) {
        i64 term = n - j + 1;
        int v_term = v_p(term, p);
        i64 u_term = term;
        for (int k = 0; k < v_term; k++) u_term /= p;
        u_term %= mod;

        int v_jd = v_p(j, p);
        i64 u_jd = j;
        for (int k = 0; k < v_jd; k++) u_jd /= p;
        u_jd %= mod;
        i64 inv_ujd = modinv(u_jd, mod);

        u_c = mmul(mmul(u_c, u_term, mod), inv_ujd, mod);
        v_c = v_c + v_term - v_jd;

        if (v_c >= a) {
            out[j] = 0;
        } else {
            i64 pp = 1;
            for (int k = 0; k < v_c; k++) pp *= p;
            out[j] = mmul(pp % mod, u_c, mod);
        }
    }
}

static i64 solve_mod_prime_power(i64 n, int p, int a) {
    i64 mod = 1;
    for (int i = 0; i < a; i++) mod *= p;

    if (n == 0) return 1 % mod;

    int J = max_j(n, p, a);
    /* J is already bounded by max_j; no need to cap at n (which overflows int) */

    i64 *fact_stirling = (i64 *)malloc((size_t)(J + 1) * sizeof(i64));
    i64 *choose = (i64 *)malloc((size_t)(J + 1) * sizeof(i64));

    forward_diff(n, J, mod, p, fact_stirling);
    binom_prefix(n, J, mod, p, a, choose);

    i64 pow2 = mpow(2, n, mod);
    i64 inv2 = modinv(2, mod);

    i64 res = 0;
    i64 cur_pow2 = pow2;
    for (int j = 0; j <= J; j++) {
        i64 term = mmul(mmul(choose[j], fact_stirling[j], mod), cur_pow2, mod);
        res = (res + term) % mod;
        cur_pow2 = mmul(cur_pow2, inv2, mod);
    }

    free(fact_stirling);
    free(choose);
    return res;
}

static i64 crt(i64 *residues, i64 *moduli, int count) {
    i64 M = 1;
    for (int i = 0; i < count; i++) M *= moduli[i];
    i128 x = 0;
    for (int i = 0; i < count; i++) {
        i64 Mi = M / moduli[i];
        i64 inv = modinv(Mi % moduli[i], moduli[i]);
        x += (i128)residues[i] * Mi * inv;
    }
    return (i64)(x % M);
}

i64 p830_native(void) {
    int primes[] = {83, 89, 97};
    int power = 3;
    i64 n = 1;
    for (int i = 0; i < 18; i++) n *= 10; /* 10^18 */

    i64 mods[3], residues[3];
    for (int i = 0; i < 3; i++) {
        mods[i] = 1;
        for (int k = 0; k < power; k++) mods[i] *= primes[i];
        residues[i] = solve_mod_prime_power(n, primes[i], power);
    }

    return crt(residues, mods, 3);
}
