#include <stdint.h>
#include <stdlib.h>

typedef int64_t i64;
typedef __int128 i128;

static const i64 MOD = 1000000007LL;
static const i64 INV2 = 500000004LL;
static const i64 INV8 = 125000001LL;

static i64 mm(i64 a, i64 b) {
    return (i64)((i128)a * b % MOD);
}

static i64 mpow(i64 a, i64 e) {
    i64 r = 1; a %= MOD; if (a < 0) a += MOD;
    while (e > 0) {
        if (e & 1) r = mm(r, a);
        a = mm(a, a);
        e >>= 1;
    }
    return r;
}

static i64 madd(i64 a, i64 b) { a += b; if (a >= MOD) a -= MOD; return a; }
static i64 msub(i64 a, i64 b) { a -= b; if (a < 0) a += MOD; return a; }

/* Returns f(n), diag_fix(n), n! mod MOD via *_out */
static void f_diag_fact(i64 n, i64 *f_out, i64 *diag_out, i64 *fact_out) {
    if (n == 0) { *f_out = 1; *diag_out = 1; *fact_out = 1; return; }

    i64 fact = 1;
    i64 h_im2 = 1, h_im1 = 0;

    i64 d0, d1, d2, d3;
    if (n <= 3) { d0 = d1 = d2 = d3 = 0; }
    else { d0 = 1; d1 = 0; d2 = 1; d3 = 4; }

    for (i64 i = 1; i <= n; i++) {
        fact = mm(fact, i % MOD);

        if (i >= 2) {
            i64 k = i - 1;
            i64 h_i = madd(mm(k % MOD, h_im1), mm(mm(k % MOD, INV2), h_im2));
            h_im2 = h_im1; h_im1 = h_i;
        }

        if (n >= 4 && i >= 4) {
            i64 k = i - 1;
            i64 term1 = mm((2 * k) % MOD, d3);
            i64 term2 = mm(mm(k % MOD, (k - 2) % MOD), d2);
            i64 term3 = mm(mm(mm(k % MOD, (k - 1) % MOD), (k - 2) % MOD), d0);
            i64 newv = msub(msub(term1, term2), mm(term3, INV2));
            d0 = d1; d1 = d2; d2 = d3; d3 = newv;
        }
    }

    *f_out = mm(fact, h_im1);

    i64 diag;
    if (n == 0) diag = 1;
    else if (n == 1) diag = 0;
    else if (n == 2) diag = 1;
    else if (n == 3) diag = 4;
    else diag = d3;
    *diag_out = diag;
    *fact_out = fact;
}

static i64 fix_axis_reflection(i64 n, i64 fact_n) {
    if (n & 1) return 0;
    return mm(fact_n, mpow(INV2, n / 2));
}

static i64 fix_rotation_90(i64 n) {
    if (n & 1) return 0;
    i64 m = n / 2;
    if (m == 0) return 1;
    if (m == 1) return 1;
    if (m == 2) return 2;

    i64 b0 = 1, b1 = 1, b2 = 2;
    for (i64 i = 2; i < m; i++) {
        i64 val = msub(madd(mm((2 * i + 1) % MOD, b2),
                            mm((2 * i % MOD) * ((i - 1) % MOD) % MOD, b0)),
                        mm(i % MOD, b1));
        b0 = b1; b1 = b2; b2 = val;
    }
    return b2;
}

static i64 fix_rotation_180(i64 n) {
    if (n == 0) return 1;

    if ((n & 1) == 0) {
        i64 m = n / 2;
        if (m == 0) return 1;
        if (m == 1) return 1;
        i64 j_prev = 1, j_curr = 1, fact = 1;
        for (i64 i = 1; i < m; i++) {
            fact = mm(fact, i % MOD);
            i64 j_next = madd(mm((4 * i + 1) % MOD, j_curr),
                              mm((4 * i) % MOD, j_prev));
            j_prev = j_curr; j_curr = j_next;
        }
        i64 fact_m = mm(fact, m % MOD);
        return mm(fact_m, j_curr);
    }

    /* odd n = 2m+1 */
    i64 m = (n - 1) / 2;
    if (m == 0) return 0;
    i64 fact = 1, t = 0, j_prev = 1, j_curr = 1;
    for (i64 i = 1; i <= m; i++) {
        fact = mm(fact, i % MOD);
        t = madd(mm((4 * i) % MOD, t), mm((2 * i) % MOD, j_prev));
        if (i < m) {
            i64 j_next = madd(mm((4 * i + 1) % MOD, j_curr),
                              mm((4 * i) % MOD, j_prev));
            j_prev = j_curr; j_curr = j_next;
        }
    }
    return mm(fact, t);
}

static i64 g(i64 n) {
    i64 f_n, diag, fact_n;
    f_diag_fact(n, &f_n, &diag, &fact_n);
    i64 axis = fix_axis_reflection(n, fact_n);
    i64 r180 = fix_rotation_180(n);
    i64 r90 = fix_rotation_90(n);

    i64 total = f_n;
    total = madd(total, r180);
    total = madd(total, mm(2, r90));
    total = madd(total, mm(2, axis));
    total = madd(total, mm(2, diag));
    return mm(total, INV8);
}

long long p741_native(void) {
    i64 n1 = 1;
    for (int i = 0; i < 7; i++) n1 *= 7;   /* 7^7 = 823543 */
    i64 n2 = 1;
    for (int i = 0; i < 8; i++) n2 *= 8;   /* 8^8 = 16777216 */
    i64 ans = madd(g(n1), g(n2));
    return (long long)ans;
}
