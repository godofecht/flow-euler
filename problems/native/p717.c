#include <stdint.h>
#include <stdlib.h>
typedef long long i64;
typedef __int128 i128;
static i64 modpow(i64 b, i64 e, i64 m) {
    i64 r = 1 % m; b %= m;
    while (e > 0) {
        if (e & 1) r = (i64)((i128)r * b % m);
        b = (i64)((i128)b * b % m);
        e >>= 1;
    }
    return r;
}
long long pe_solve(void) {
    int N = 10000000;
    char *comp = calloc((size_t)N + 1, 1);
    for (int i = 2; (i64)i * i <= N; i++)
        if (!comp[i])
            for (int j = i * i; j <= N; j += i) comp[j] = 1;
    i64 total = 0;
    for (int p = 3; p <= N; p++) {
        if (comp[p]) continue;
        i64 x = modpow(2, p, p - 1);
        i64 r = modpow(2, x, p);
        i64 n = r * (p - 1) / (2 * p) + 1;
        i64 modpp = (i64)p * p;
        i64 pp = modpow(2, p - 1, modpp);
        i64 t = (i64)(((-((i128)r * (p - 1) % modpp) * pp) % modpp + 1) % modpp);
        if (t < 0) t += modpp;
        i64 fp = t / p;
        total += (fp + 2 * n) % p;
    }
    free(comp);
    return total;
}
