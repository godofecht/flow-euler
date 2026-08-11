#include <stdint.h>
#include <string.h>
#include <stdlib.h>
typedef long long i64;

/*
 * Project Euler 768: Chandelier.
 *
 * Count arrangements of 20 identical candles in 360 distinct sockets such that
 * the vector sum of the chosen 360th roots of unity is 0.
 *
 * Factorisation 360 = 5 * 8 * 9 and cyclotomic structure reduce the balance
 * constraint to 12 identical independent blocks.  For each block we build a
 * generating function S(x) by summing (P_delta(x))^3 over all differences delta
 * of two pentagon-choices, then the answer is [x^20] S(x)^12.
 *
 * Exact integer arithmetic is recovered via CRT over three 64-bit primes.
 */

#define M 20
#define N (M + 1)  /* polynomial size: degrees 0..20 */

static const i64 P1 = 998244353LL, P2 = 1004535809LL, P3 = 469762049LL;

static i64 madd(i64 a, i64 b, i64 m) { a += b; if (a >= m) a -= m; return a; }
static i64 msub(i64 a, i64 b, i64 m) { a -= b; if (a < 0) a += m; return a; }
static i64 mmul(i64 a, i64 b, i64 m) { return (i64)((__int128)a * b % m); }
static i64 mpow(i64 a, i64 e, i64 m) {
    i64 r = 1 % m; a %= m;
    while (e > 0) { if (e & 1) r = mmul(r, a, m); a = mmul(a, a, m); e >>= 1; }
    return r;
}
static i64 minv(i64 a, i64 m) { return mpow(a, m - 2, m); }

/* res = a * b  truncated to degree M (mod m) */
static void polymul(i64 res[N], const i64 a[N], const i64 b[N], i64 mod) {
    i64 tmp[N];
    memset(tmp, 0, sizeof(tmp));
    for (int i = 0; i < N; i++) {
        if (a[i] == 0) continue;
        for (int j = 0; j < N - i; j++) {
            if (b[j] == 0) continue;
            tmp[i + j] = madd(tmp[i + j], mmul(a[i], b[j], mod), mod);
        }
    }
    memcpy(res, tmp, sizeof(tmp));
}

/* res = base^exp  truncated to degree M (mod m) */
static void polypow(i64 res[N], const i64 base[N], int exp, i64 mod) {
    memset(res, 0, sizeof(i64) * N);
    res[0] = 1 % mod;
    i64 cur[N];
    memcpy(cur, base, sizeof(cur));
    while (exp > 0) {
        if (exp & 1) polymul(res, res, cur, mod);
        exp >>= 1;
        if (exp) polymul(cur, cur, cur, mod);
    }
}

typedef struct { int key; int w; } PairEntry;

static int cmp_pair(const void *a, const void *b) {
    const PairEntry *pa = a, *pb = b;
    if (pa->key != pb->key) return pa->key < pb->key ? -1 : 1;
    return 0;
}

/* Solve f(360,20) modulo one prime. */
static i64 solve_mod(i64 mod) {
    /* Build the 32 pentagon-choices (subsets of {1,y,y^2,y^3,y^4}). */
    i64 pat_coeff[32][4];
    int  pat_cnt[32];
    for (int mask = 0; mask < 32; mask++) {
        i64 coeff[4] = {0, 0, 0, 0};
        int cnt = 0;
        for (int v = 0; v < 5; v++) {
            if (mask & (1 << v)) {
                cnt++;
                if (v < 4) coeff[v] += 1;
                else { for (int k = 0; k < 4; k++) coeff[k] -= 1; } /* y^4 = -(1+y+y^2+y^3) */
            }
        }
        for (int k = 0; k < 4; k++) pat_coeff[mask][k] = coeff[k];
        pat_cnt[mask] = cnt;
    }

    /* Collect all 32*32 ordered pairs, keyed by their 4-tuple difference. */
    PairEntry pairs[1024];
    int npairs = 0;
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 32; j++) {
            int d[4];
            for (int k = 0; k < 4; k++)
                d[k] = (int)(pat_coeff[i][k] - pat_coeff[j][k]) + 10;  /* shift to [0,20] */
            int key = d[0] + 21 * d[1] + 441 * d[2] + 9261 * d[3];
            pairs[npairs].key = key;
            pairs[npairs].w   = pat_cnt[i] + pat_cnt[j];
            npairs++;
        }
    }
    qsort(pairs, npairs, sizeof(PairEntry), cmp_pair);

    /* Group by delta, cube each P_delta, accumulate into S. */
    i64 S[N];
    memset(S, 0, sizeof(S));

    int idx = 0;
    while (idx < npairs) {
        int key = pairs[idx].key;
        i64 P[N];
        memset(P, 0, sizeof(P));
        while (idx < npairs && pairs[idx].key == key) {
            int w = pairs[idx].w;
            if (w <= M) P[w] = madd(P[w], 1, mod);
            idx++;
        }
        i64 P2[N], P3[N];
        polymul(P2, P, P, mod);
        polymul(P3, P2, P, mod);
        for (int i = 0; i < N; i++) S[i] = madd(S[i], P3[i], mod);
    }

    /* 12 independent blocks: answer = [x^20] S(x)^12. */
    i64 result[N];
    polypow(result, S, 12, mod);
    return result[M];
}

/* CRT over three primes to recover the exact integer. */
static i64 crt3(i64 a1, i64 a2, i64 a3) {
    i64 inv12 = minv(P1 % P2, P2);
    __int128 P12 = (__int128)P1 * P2;
    i64 t1 = msub(a2, a1 % P2, P2);
    t1 = mmul(t1, inv12, P2);
    __int128 x2 = (__int128)a1 + (__int128)P1 * t1;
    i64 inv123 = minv((i64)(P12 % P3), P3);
    i64 t2 = msub(a3, (i64)(x2 % P3), P3);
    t2 = mmul(t2, inv123, P3);
    __int128 x = x2 + P12 * t2;
    return (i64)x;
}

long long p768_native(void) {
    i64 a1 = solve_mod(P1);
    i64 a2 = solve_mod(P2);
    i64 a3 = solve_mod(P3);
    return crt3(a1, a2, a3);
}
