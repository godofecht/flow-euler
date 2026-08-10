#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef __int128 i128;
typedef long long i64;

static i64 rng = 123456789;
static i64 rnd(i64 modv) {
    rng = rng * 6364136223846793005LL + 1;
    i64 v = rng;
    if (v < 0) v = -v;
    return v % modv;
}

static i128 igcd(i128 a, i128 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { i128 t = a % b; a = b; b = t; }
    return a;
}

static i128 modmul(i128 a, i128 b, i128 n) { return a * b % n; }

static i128 modpow(i128 base, i128 exp, i128 n) {
    i128 r = 1, b = base % n, e = exp;
    while (e) {
        if (e & 1) r = modmul(r, b, n);
        b = modmul(b, b, n);
        e >>= 1;
    }
    return r;
}

static int is_prime(i128 n) {
    if (n < 2) return 0;
    i64 sp[] = {2,3,5,7,11,13,17,19,23,29,0};
    for (int i=0;sp[i];i++) if (n % sp[i] == 0) return n == sp[i];
    i128 d = n - 1; i64 s = 0;
    while ((d & 1) == 0) { s++; d >>= 1; }
    i64 as[] = {2,325,9375,28178,450775,9780504,1795265022,0};
    for (int i=0;as[i];i++) {
        i128 a = as[i] % n; if (!a) continue;
        i128 x = modpow(a, d, n);
        if (x == 1 || x == n - 1) continue;
        int ok = 0;
        for (i64 j=1;j<s;j++) {
            x = modmul(x, x, n);
            if (x == n - 1) { ok = 1; break; }
        }
        if (!ok) return 0;
    }
    return 1;
}

static i128 pollard(i128 n) {
    if ((n & 1) == 0) return 2;
    if (n % 3 == 0) return 3;
    for (;;) {
        i128 x = 2 + rnd(1000000007) % n;
        i128 y = x;
        i128 c = 1 + rnd(1000000007) % n;
        i128 d = 1;
        while (d == 1) {
            x = (modmul(x, x, n) + c) % n;
            y = (modmul(y, y, n) + c) % n;
            y = (modmul(y, y, n) + c) % n;
            i128 diff = x > y ? x - y : y - x;
            d = igcd(diff, n);
        }
        if (d != n) return d;
    }
}

static void factor_rec(i128 n, i128 *ps, i64 *es, int *k) {
    if (n == 1) return;
    if (is_prime(n)) {
        for (int i=0;i<*k;i++) if (ps[i]==n) { es[i]++; return; }
        ps[*k]=n; es[*k]=1; (*k)++;
        return;
    }
    i128 d = pollard(n);
    factor_rec(d, ps, es, k);
    factor_rec(n / d, ps, es, k);
}

static void next_triangle_jump(i128 a_tri, i128 *out_k, i128 *out_a) {
    i128 M = 2 * a_tri;
    i128 ps[64]; i64 es[64]; int fk = 0;
    factor_rec(M, ps, es, &fk);
    i128 *divs = malloc(50000 * sizeof(i128));
    int nd = 1; divs[0] = 1;
    for (int i=0;i<fk;i++) {
        int old = nd; i128 pe = 1;
        for (i64 e=0;e<es[i];e++) {
            pe *= ps[i];
            for (int j=0;j<old;j++) divs[nd++] = divs[j] * pe;
        }
    }
    for (int i=0;i<nd;i++) for (int j=i+1;j<nd;j++) if (divs[j]<divs[i]) {
        i128 t=divs[i]; divs[i]=divs[j]; divs[j]=t;
    }
    for (int i=0;i<nd;i++) {
        i128 x = divs[i], y = M / x;
        if ((x - y) & 1) {
            i128 k = (x - y - 1) / 2;
            if (k > 0) {
                *out_k = k;
                *out_a = a_tri + k * (k + 1) / 2;
                free(divs);
                return;
            }
        }
    }
    free(divs);
    *out_k = -1; *out_a = -1;
}

long long pe_solve(void) {
    i128 n = 0, a = 3;
    for (int i=0;i<69;i++) {
        i128 step, na;
        next_triangle_jump(a, &step, &na);
        n += step; a = na;
    }
    return (long long)n;
}
