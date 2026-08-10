#include <stdlib.h>
#include <string.h>
#include <gmp.h>

static long comb_l(long n, long k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    long r = 1;
    for (long i = 1; i <= k; i++) r = r * (n - k + i) / i;
    return r;
}

double pe687_answer(void) {
    const int RANKS = 13, N = 52, BASE = 24, MAXD = 40;
    mpz_t Qpow[14][40];
    int deg[14];
    for (int m = 0; m <= RANKS; m++)
        for (int i = 0; i < MAXD; i++) mpz_init_set_ui(Qpow[m][i], 0);
    mpz_set_ui(Qpow[0][0], 1);
    deg[0] = 0;
    long Q[4] = {1, -12, 36, -24};
    for (int m = 1; m <= RANKS; m++) {
        int pd = deg[m - 1];
        deg[m] = pd + 3;
        for (int i = 0; i <= pd; i++) {
            for (int j = 0; j < 4; j++) {
                mpz_t tmp; mpz_init(tmp);
                mpz_mul_si(tmp, Qpow[m - 1][i], Q[j]);
                mpz_add(Qpow[m][i + j], Qpow[m][i + j], tmp);
                mpz_clear(tmp);
            }
        }
    }

    mpz_t fact[53], denom, num, nval[14], z, xk, good, tmp, total, scale, q;
    for (int i = 0; i <= N; i++) mpz_init(fact[i]);
    mpz_set_ui(fact[0], 1);
    for (int i = 1; i <= N; i++) mpz_mul_ui(fact[i], fact[i - 1], (unsigned)i);
    mpz_init_set_ui(denom, 1);
    for (int i = 0; i < RANKS; i++) mpz_mul_ui(denom, denom, BASE);
    mpz_inits(num, z, xk, good, tmp, total, scale, q, NULL);
    for (int m = 0; m <= RANKS; m++) mpz_init(nval[m]);

    for (int m = 0; m <= RANKS; m++) {
        mpz_set_ui(num, 0);
        for (int B = 0; B <= deg[m]; B++) {
            mpz_mul(tmp, fact[N - B], Qpow[m][B]);
            mpz_add(num, num, tmp);
        }
        mpz_divexact(nval[m], num, denom);
    }
    mpz_set(total, nval[0]);

    mpz_set_ui(good, 0);
    int primes[] = {2, 3, 5, 7, 11, 13};
    for (int ki = 0; ki < 6; ki++) {
        int k = primes[ki];
        mpz_set_ui(z, 0);
        for (int m = k; m <= RANKS; m++) {
            long c = comb_l(RANKS - k, m - k);
            if ((m - k) & 1) c = -c;
            mpz_mul_si(tmp, nval[m], c);
            mpz_add(z, z, tmp);
        }
        mpz_mul_ui(xk, z, (unsigned long)comb_l(RANKS, k));
        mpz_add(good, good, xk);
    }

    mpz_ui_pow_ui(scale, 10, 11);
    mpz_mul(q, good, scale);
    mpz_fdiv_q_ui(tmp, total, 2);
    mpz_add(q, q, tmp);
    mpz_fdiv_q(q, q, total);
    mpz_add_ui(q, q, 5);
    mpz_fdiv_q_ui(q, q, 10);
    double ans = mpz_get_d(q) / 1e10;

    for (int m = 0; m <= RANKS; m++)
        for (int i = 0; i < MAXD; i++) mpz_clear(Qpow[m][i]);
    for (int i = 0; i <= N; i++) mpz_clear(fact[i]);
    for (int m = 0; m <= RANKS; m++) mpz_clear(nval[m]);
    mpz_clears(denom, num, z, xk, good, tmp, total, scale, q, NULL);
    return ans;
}
