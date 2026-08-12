#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/* Project Euler 824 - Chess Sliders on a cylindrical board.
 *
 * L(N, K) mod (10^7+19)^2 where N=10^9, K=10^15.
 *
 * Uses Wilson quotient, harmonic numbers mod p, unit factorials mod p^2.
 * __int128 is used for modular multiplication mod p^2.
 */

static const uint64_t P = 10000019ULL;
static const uint64_t MOD = 10000019ULL * 10000019ULL;  /* p^2 */

static uint32_t *invp;
static uint64_t *fac;
static uint32_t *H;
static uint64_t w;  /* Wilson quotient */

static inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m) {
    return (uint64_t)((__uint128_t)a * b % m);
}

static uint64_t inv_mod_p2(uint64_t a) {
    a %= MOD;
    uint64_t r = a % P;
    uint64_t x = invp[r];
    uint64_t ax = mulmod(a, x, MOD);
    uint64_t term = (2 + MOD - ax) % MOD;
    return mulmod(x, term, MOD);
}

static uint64_t Fpow(uint64_t u) {
    uint64_t um = u % P;
    uint64_t a1 = um * P;               /* < P^2 = MOD */
    uint64_t a2 = mulmod(a1, w, MOD);
    uint64_t tmp = (1 + MOD - a2) % MOD;
    if (u & 1) tmp = (MOD - tmp) % MOD;
    return tmp;
}

static uint64_t unit_factorial(uint64_t n) {
    uint64_t res = 1;
    while (n) {
        uint64_t u = n / P;
        uint64_t v = n % P;
        res = mulmod(res, fac[v], MOD);
        uint64_t corr = (u % P) * H[v] % P;
        res = mulmod(res, 1 + corr * P, MOD);
        res = mulmod(res, Fpow(u), MOD);
        n = u;
    }
    return res;
}

static uint64_t vp_fact(uint64_t n) {
    uint64_t q = n / P;
    return q + q / P;
}

static uint64_t binom_mod_p2(uint64_t n, uint64_t k) {
    if (k > n) return 0;
    uint64_t nk = n - k;
    uint64_t e = vp_fact(n) - vp_fact(k) - vp_fact(nk);
    if (e >= 2) return 0;
    uint64_t un = unit_factorial(n);
    uint64_t uk = unit_factorial(k);
    uint64_t unk = unit_factorial(nk);
    uint64_t val = mulmod(un, inv_mod_p2(uk), MOD);
    val = mulmod(val, inv_mod_p2(unk), MOD);
    if (e == 1) val = mulmod(val, P, MOD);
    return val;
}

static uint64_t coeff_alpha_power(uint64_t M, uint64_t d) {
    if (d == 0) return 1;
    uint64_t B = binom_mod_p2(M - d - 1, d - 1);
    uint64_t val = mulmod(M % MOD, inv_mod_p2(d % MOD), MOD);
    return mulmod(val, B, MOD);
}

long long p824_native(void) {
    invp = calloc(P, sizeof(uint32_t));
    H    = calloc(P, sizeof(uint32_t));
    fac  = malloc(P * sizeof(uint64_t));
    if (!invp || !fac || !H) {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }

    /* Inverses mod p */
    invp[1] = 1;
    for (uint64_t i = 2; i < P; i++) {
        invp[i] = (uint32_t)((P - (P / i) * invp[P % i] % P) % P);
    }

    /* Factorial mod p^2 for 0..p-1 */
    uint64_t f = 1;
    fac[0] = 1;
    for (uint64_t i = 1; i < P; i++) {
        f = mulmod(f, i, MOD);
        fac[i] = f;
    }

    /* Harmonic numbers H[v] = sum_{i<=v} 1/i mod p */
    uint32_t h = 0;
    for (uint64_t i = 1; i < P; i++) {
        h = (uint32_t)((h + invp[i]) % P);
        H[i] = h;
    }

    /* Wilson quotient: (p-1)! = -1 + p*w (mod p^2) */
    uint64_t F = fac[P - 1];
    w = ((F + 1) / P) % P;

    /* Main computation */
    uint64_t N = 1000000000ULL;        /* 10^9 */
    uint64_t K = 1000000000000000ULL;  /* 10^15 */
    uint64_t t_max = K / N;            /* 10^6 */

    uint64_t comb = 1;
    uint64_t ans = 0;
    uint64_t d = K;
    uint64_t M = N * N;

    for (uint64_t t = 0; t <= t_max; t++) {
        if (t) {
            comb = mulmod(comb, (N - t + 1) % MOD, MOD);
            comb = mulmod(comb, inv_mod_p2(t), MOD);
        }
        ans = (ans + mulmod(comb, coeff_alpha_power(M, d), MOD)) % MOD;
        d -= N;
        M -= 2 * N;
    }

    free(invp);
    free(fac);
    free(H);

    return (long long)ans;
}
