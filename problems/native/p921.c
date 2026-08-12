#include <stdint.h>

typedef long long i64;
typedef unsigned long long u64;

/* Project Euler 921 - Golden Recurrence
   S(1618034) mod 398874989
*/

#define MOD 398874989LL
#define TARGET_M 1618034LL

static i64 powmod(i64 a, i64 b, i64 m) {
    i64 r = 1; a %= m; if (a < 0) a += m;
    while (b > 0) {
        if (b & 1) r = (__int128)r * a % m;
        a = (__int128)a * a % m;
        b >>= 1;
    }
    return r;
}

static i64 legendre(i64 a, i64 p) {
    return powmod(a, (p - 1) / 2, p);
}

static i64 sqrt_mod_prime(i64 n, i64 p) {
    n %= p; if (n < 0) n += p;
    if (n == 0) return 0;
    if (p == 2) return n;
    if (legendre(n, p) != 1) return -1;
    if (p % 4 == 3) return powmod(n, (p + 1) / 4, p);
    i64 q = p - 1; int s = 0;
    while (q % 2 == 0) { q /= 2; s++; }
    i64 z = 2;
    while (legendre(z, p) != p - 1) z++;
    i64 c = powmod(z, q, p);
    i64 x = powmod(n, (q + 1) / 2, p);
    i64 t = powmod(n, q, p);
    int m = s;
    while (t != 1) {
        i64 t2 = t; int i = 1;
        for (; i < m; i++) {
            t2 = (__int128)t2 * t2 % p;
            if (t2 == 1) break;
        }
        i64 b = powmod(c, 1LL << (m - i - 1), p);
        x = (__int128)x * b % p;
        i64 bb = (__int128)b * b % p;
        t = (__int128)t * bb % p;
        c = bb;
        m = i;
    }
    return x;
}

long long p921_native(void) {
    i64 p = MOD;
    i64 exp_mod = p - 1;
    i64 inv2 = powmod(2, p - 2, p);

    i64 r = sqrt_mod_prime(5, p);
    i64 inv_r = powmod(r, p - 2, p);

    i64 phi = ((__int128)(1 + r) * inv2) % p;
    i64 g = (__int128)phi * phi % p;
    g = (__int128)g * phi % p;

    /* Check orientation: F_3 should be 2 */
    i64 u = g;
    i64 uinv = powmod(u, p - 2, p);
    i64 f3 = (__int128)((u + uinv) % p) * inv_r % p;
    if (f3 != 2) {
        r = (p - r) % p;
        inv_r = powmod(r, p - 2, p);
        phi = ((__int128)(1 + r) * inv2) % p;
        g = (__int128)phi * phi % p;
        g = (__int128)g * phi % p;
        u = g;
        uinv = powmod(u, p - 2, p);
    }

    i64 inv32 = powmod(32, p - 2, p);

    /* Process i=2: E_1 = E_2 = 5 */
    i64 e_im2 = 5 % exp_mod;
    i64 e_im1 = 5 % exp_mod;

    i64 total = 0;

    /* i=2 */
    i64 e = e_im1;
    u = powmod(g, e, p);
    uinv = powmod(u, p - 2, p);
    i64 f = (__int128)((u + uinv) % p) * inv_r % p;
    i64 l = (u - uinv) % p; if (l < 0) l += p;
    i64 f2 = (__int128)f * f % p;
    i64 f4 = (__int128)f2 * f2 % p;
    i64 f5 = (__int128)f4 * f % p;
    i64 l2 = (__int128)l * l % p;
    i64 l4 = (__int128)l2 * l2 % p;
    i64 l5 = (__int128)l4 * l % p;
    i64 term = (f5 + l5) % p;
    term = (__int128)term * inv32 % p;
    total = term;

    /* i=3..m */
    for (i64 i = 3; i <= TARGET_M; i++) {
        e = (__int128)e_im1 * e_im2 % exp_mod;
        e_im2 = e_im1;
        e_im1 = e;

        u = powmod(g, e, p);
        uinv = powmod(u, p - 2, p);

        f = (__int128)((u + uinv) % p) * inv_r % p;
        l = (u - uinv) % p; if (l < 0) l += p;

        f2 = (__int128)f * f % p;
        f4 = (__int128)f2 * f2 % p;
        f5 = (__int128)f4 * f % p;
        l2 = (__int128)l * l % p;
        l4 = (__int128)l2 * l2 % p;
        l5 = (__int128)l4 * l % p;

        term = (f5 + l5) % p;
        term = (__int128)term * inv32 % p;

        total += term;
        if (total >= p) total -= p;
    }

    return total;
}
