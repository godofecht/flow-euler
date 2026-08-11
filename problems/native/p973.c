#include <stdint.h>
#include <stdlib.h>

#define MOD 1000000007LL
#define INV2 500000004LL
#define N 10000

static long long mod(long long x) {
    x %= MOD;
    if (x < 0) x += MOD;
    return x;
}

static long long powmod(long long base, long long exp, long long m) {
    long long result = 1 % m;
    base %= m;
    if (base < 0) base += m;
    while (exp > 0) {
        if (exp & 1) result = result * base % m;
        base = base * base % m;
        exp >>= 1;
    }
    return result;
}

static long long S_k_value(int n, int k) {
    if (k == 0) {
        if (n == 0) return 1;
        return mod(-powmod(MOD - 2, n - 1, MOD));
    }
    long long B = 1LL << k;
    long long period = B << 1;
    long long per_mask = period - 1;
    long long *S = (long long *)calloc((size_t)(n + 1), sizeof(long long));
    long long *P = (long long *)calloc((size_t)(n + 1), sizeof(long long));
    S[0] = 1;
    P[0] = 1;
    for (int nn = 1; nn <= n; nn++) {
        long long val = 0;
        long long s = 1;
        while (s <= nn) {
            long long p = s & per_mask;
            long long sign, e;
            if (p < B) {
                sign = 1;
                e = s + (B - p) - 1;
            } else {
                sign = -1;
                e = s + (period - p) - 1;
            }
            if (e > nn) e = nn;
            long long t_lo = nn - e;
            long long t_hi = nn - s;
            long long sum_range;
            if (t_lo == 0) {
                sum_range = P[t_hi];
            } else {
                sum_range = P[t_hi] - P[t_lo - 1];
            }
            val += sign * sum_range;
            s = e + 1;
        }
        S[nn] = mod(val);
        P[nn] = mod(P[nn - 1] + S[nn]);
    }
    long long result = S[n];
    free(S);
    free(P);
    return result;
}

static long long X(int n) {
    if (n <= 0) return 0;
    long long P = powmod(2, (long long)(n - 1), MOD);
    long long total = 0;
    for (int k = 0; k < 14; k++) {
        long long Sk = S_k_value(n, k);
        long long Ak = mod(((P - Sk) % MOD) * INV2 % MOD);
        total = (total + powmod(2, (long long)k, MOD) * Ak) % MOD;
    }
    total = mod(total - (long long)(n & 1));
    return total;
}

long long pe973_answer(void) {
    return X(N);
}
