#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef long long i64;
typedef unsigned long long u64;
typedef __uint128_t u128;

#define MOD 993353399LL

/* ---- 64-bit modular arithmetic (for Miller-Rabin / Pollard-Rho) ---- */

static u64 mulmod(u64 a, u64 b, u64 m) {
    return (u64)((u128)a * b % m);
}

static u64 powmod(u64 a, u64 e, u64 m) {
    u64 r = 1 % m;
    a %= m;
    while (e) {
        if (e & 1) r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}

/* ---- deterministic Miller-Rabin for 64-bit ---- */

static int is_prime_64(u64 n) {
    if (n < 2) return 0;
    static const u64 sp[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++) {
        if (n == sp[i]) return 1;
        if (n % sp[i] == 0) return 0;
    }
    u64 d = n - 1;
    int s = 0;
    while ((d & 1) == 0) { d >>= 1; s++; }
    static const u64 bases[] = {2,325,9375,28178,450775,9780504,1795265022ULL};
    for (int i = 0; i < 7; i++) {
        u64 a = bases[i];
        if (a % n == 0) continue;
        u64 x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        int comp = 1;
        for (int j = 0; j < s - 1; j++) {
            x = mulmod(x, x, n);
            if (x == n - 1) { comp = 0; break; }
        }
        if (comp) return 0;
    }
    return 1;
}

/* ---- Pollard-Rho ---- */

static u64 rng_state = 0x9E3779B97F4A7C15ULL;
static u64 rand64(void) {
    u64 x = rng_state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    rng_state = x;
    return x * 2685821657736338717ULL;
}

static u64 gcd_u64(u64 a, u64 b) {
    while (b) { u64 t = a % b; a = b; b = t; }
    return a;
}

static u64 pollard_rho(u64 n) {
    if ((n & 1) == 0) return 2;
    if (n % 3 == 0) return 3;
    while (1) {
        u64 c = rand64() % (n - 1) + 1;
        u64 x = rand64() % (n - 2) + 2;
        u64 y = x;
        u64 d = 1;
        while (d == 1) {
            x = (mulmod(x, x, n) + c) % n;
            y = (mulmod(y, y, n) + c) % n;
            y = (mulmod(y, y, n) + c) % n;
            u64 diff = x > y ? x - y : y - x;
            d = gcd_u64(diff, n);
        }
        if (d != n) return d;
    }
}

/* ---- factorization ---- */

#define MAXF 64
typedef struct { u64 prime; int exp; } Factor;

static void factorize(u64 n, Factor *out, int *outn) {
    *outn = 0;
    if (n <= 1) return;
    /* trial division by integers up to 10000 */
    for (u64 p = 2; p * p <= n && p <= 10000; p++) {
        if (n % p == 0) {
            int e = 0;
            while (n % p == 0) { n /= p; e++; }
            out[*outn].prime = p; out[*outn].exp = e; (*outn)++;
        }
    }
    if (n == 1) return;
    if (is_prime_64(n)) {
        out[*outn].prime = n; out[*outn].exp = 1; (*outn)++;
        return;
    }
    /* Pollard-Rho with explicit stack */
    u64 stack[64]; int sp = 0;
    stack[sp++] = n;
    while (sp) {
        u64 m = stack[--sp];
        if (m == 1) continue;
        if (is_prime_64(m)) {
            int found = 0;
            for (int i = 0; i < *outn; i++) {
                if (out[i].prime == m) { out[i].exp++; found = 1; break; }
            }
            if (!found) { out[*outn].prime = m; out[*outn].exp = 1; (*outn)++; }
            continue;
        }
        u64 d = pollard_rho(m);
        stack[sp++] = d;
        stack[sp++] = m / d;
    }
}

/* ---- modular arithmetic for the answer modulus (fits in 32 bits) ---- */

static i64 powmod_small(i64 a, long long e, i64 m) {
    /* e >= 0 always for our usage */
    i64 r = 1 % m;
    a %= m;
    if (a < 0) a += m;
    while (e > 0) {
        if (e & 1) r = r * a % m;
        a = a * a % m;
        e >>= 1;
    }
    return r;
}

/* ---- g(q^e) mod ---- */
/*
   g(q^e) = (q-1)^3 * sum_{t=1..e} t^2 * q^{3e-t-2}
            + q^{2e-2} * (e*(q-1)+q)^2
*/
static i64 g_prime_power_mod(u64 q, int e, i64 mod) {
    i64 qm = (i64)(q % (u64)mod);
    i64 s = 0;
    for (int t = 1; t <= e; t++) {
        long long exp = 3LL * e - t - 2;
        i64 term = powmod_small(qm, exp, mod);
        s = (s + (i64)(t * t) % mod * term) % mod;
    }
    i64 qm1 = (i64)((q - 1) % (u64)mod);
    i64 term1 = powmod_small(qm1, 3, mod) * s % mod;
    i64 inner = (i64)(((u64)e * (q - 1) + q) % (u64)mod);
    i64 term2 = powmod_small(qm, 2 * e - 2, mod) * powmod_small(inner, 2, mod) % mod;
    return (term1 + term2) % mod;
}

static i64 g_from_factorization(Factor *facs, int nf, i64 mod) {
    i64 g = 1 % mod;
    for (int i = 0; i < nf; i++)
        g = g * g_prime_power_mod(facs[i].prime, facs[i].exp, mod) % mod;
    return g;
}

static i64 f_of_prime(u64 p, i64 mod) {
    u64 m = p - 1;
    Factor facs[MAXF]; int nf;
    factorize(m, facs, &nf);
    i64 g = g_from_factorization(facs, nf, mod);
    i64 mm = (i64)(m % (u64)mod);
    return (mm * mm + g) % mod;
}

/* ---- small prime sieve up to limit ---- */

static int sieve_small(int limit, int *primes) {
    static char comp[200005];
    memset(comp, 0, (size_t)limit + 1);
    int pc = 0;
    for (int p = 2; p <= limit; p++) {
        if (!comp[p]) {
            primes[pc++] = p;
            for (long long j = (long long)p * p; j <= limit; j += p)
                comp[j] = 1;
        }
    }
    return pc;
}

/* ---- main solver ---- */

long long p801_native(void) {
    i64 mod = MOD;
    u64 lo = 10000000000000000ULL;   /* 10^16 */
    u64 hi = lo + 1000000ULL;         /* 10^16 + 10^6 */
    long long length = (long long)(hi - lo + 1);

    char *is_comp = (char *)calloc((size_t)length, 1);
    if (!is_comp) { fprintf(stderr, "oom\n"); exit(1); }

    int primes[20000];
    int pc = sieve_small(200000, primes);

    for (int i = 0; i < pc; i++) {
        u64 q = (u64)primes[i];
        u64 offset = (q - (lo % q)) % q;
        for (long long j = (long long)offset; j < length; j += (long long)q)
            is_comp[j] = 1;
        if (lo <= q && q <= hi)
            is_comp[q - lo] = 0;
    }

    i64 total = 0;
    for (long long j = 0; j < length; j++) {
        if (is_comp[j]) continue;
        u64 n = lo + (u64)j;
        if (is_prime_64(n))
            total = (total + f_of_prime(n, mod)) % mod;
    }

    free(is_comp);
    return total;
}
