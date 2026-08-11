#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef long long i64;

#define MOD 1000000007LL
#define M 200

static i64 gcd_ll(i64 a, i64 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { i64 t = a % b; a = b; b = t; }
    return a;
}

static i64 mod_pow(i64 a, i64 e, i64 mod) {
    i64 r = 1 % mod; a %= mod; if (a < 0) a += mod;
    while (e > 0) {
        if (e & 1) r = (i64)((__int128)r * a % mod);
        a = (i64)((__int128)a * a % mod);
        e >>= 1;
    }
    return r;
}

static i64 mod_inv(i64 a, i64 mod) { return mod_pow(a, mod - 2, mod); }

/* Sieve primes up to 10000 (sqrt(80e6) < 9000). */
static int primes[2000];
static int nprimes = 0;

static void init_primes(void) {
    int limit = 10000;
    char *sieve = calloc(limit + 1, 1);
    for (int i = 2; i <= limit; i++) {
        if (!sieve[i]) {
            primes[nprimes++] = i;
            for (int j = i * i; j <= limit; j += i) sieve[j] = 1;
        }
    }
    free(sieve);
}

/* Factorize n into arrays of prime and exponent. Returns count. */
static int factorize(i64 n, i64 *pp, int *ee) {
    int cnt = 0;
    i64 x = n;
    for (int i = 0; i < nprimes; i++) {
        i64 p = primes[i];
        if (p * p > x) break;
        if (x % p == 0) {
            int e = 0;
            while (x % p == 0) { x /= p; e++; }
            pp[cnt] = p; ee[cnt] = e; cnt++;
        }
    }
    if (x > 1) {
        pp[cnt] = x; ee[cnt] = 1; cnt++;
    }
    return cnt;
}

static i64 euler_phi(i64 n) {
    if (n == 1) return 1;
    i64 pp[64]; int ee[64];
    int c = factorize(n, pp, ee);
    i64 phi = n;
    for (int i = 0; i < c; i++) phi = (phi / pp[i]) * (pp[i] - 1);
    return phi;
}

/* Multiplicative order of 10 mod m. Returns 0 if gcd(10,m) != 1. */
static i64 mult_order_10(i64 m) {
    if (m == 1) return 1;
    if (gcd_ll(m, 10) != 1) return 0;
    i64 phi = euler_phi(m);
    i64 pp[64]; int ee[64];
    int c = factorize(phi, pp, ee);
    i64 k = phi;
    for (int i = 0; i < c; i++) {
        i64 p = pp[i];
        while (k % p == 0 && mod_pow(10, k / p, m) == 1) {
            k /= p;
        }
    }
    return k;
}

/* Derive bounds on digit length k for leading digit d.
   k_low = -1 means no solution. k_high = -1 means unbounded. */
static void digit_k_bounds(i64 a, i64 b, i64 d, i64 *k_low, i64 *k_high) {
    i64 e = 0, t = a;
    while (t < b) { t *= 10; e++; }
    *k_low = e + 1;

    i64 denom = a * (d + 1) - 10 * b;
    if (denom <= 0) {
        *k_high = -1;
        return;
    }
    i64 max_pow10 = (b * d - 1) / denom;
    if (max_pow10 <= 0) {
        *k_low = -1;
        *k_high = -1;
        return;
    }
    i64 e_high = 0, pow10 = 1;
    while (pow10 * 10 <= max_pow10) {
        pow10 *= 10;
        e_high++;
    }
    *k_high = e_high + 1;
}

/* For reduced a/b, find minimal (k,d,D). Returns 1 on success. */
static int find_N_params(i64 a, i64 b, i64 *out_k, i64 *out_d, i64 *out_D) {
    if (a == b) {
        *out_k = 1; *out_d = 1; *out_D = 9 * b;
        return 1;
    }
    i64 D = 10 * b - a;
    if (D <= 0) return 0;

    i64 best_k = 0, best_d = 0; int have = 0;
    for (i64 d = 1; d <= 9; d++) {
        i64 k_low, k_high;
        digit_k_bounds(a, b, d, &k_low, &k_high);
        if (k_low == -1) continue;
        if (k_low < 2) k_low = 2;

        i64 m = D / gcd_ll(D, d * b);
        i64 ord10 = mult_order_10(m);
        if (ord10 == 0) continue;

        i64 k = ((k_low + ord10 - 1) / ord10) * ord10;
        if (k_high != -1 && k > k_high) continue;

        if (!have || k < best_k || (k == best_k && d < best_d)) {
            best_k = k; best_d = d; have = 1;
        }
    }
    if (!have) return 0;
    *out_k = best_k; *out_d = best_d; *out_D = D;
    return 1;
}

static i64 N_mod(i64 a, i64 b) {
    i64 g = gcd_ll(a, b);
    a /= g; b /= g;
    if (a == b) return 1 % MOD;

    i64 k, d, D;
    if (!find_N_params(a, b, &k, &d, &D)) return 0;

    i64 invD = mod_inv(D % MOD, MOD);
    i64 ten_k = mod_pow(10, k, MOD);
    i64 term = (ten_k - 1 + MOD) % MOD;
    i64 res = (d % MOD) * (b % MOD) % MOD;
    res = res * term % MOD;
    res = res * invD % MOD;
    return res;
}

long long p805_native(void) {
    init_primes();
    i64 total = 0;
    for (int u = 1; u <= M; u++) {
        i64 a = (i64)u * u * u;
        for (int v = 1; v <= M; v++) {
            if (gcd_ll(u, v) != 1) continue;
            i64 b = (i64)v * v * v;
            total = (total + N_mod(a, b)) % MOD;
        }
    }
    return total;
}
