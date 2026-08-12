#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;
typedef uint64_t u64;
typedef uint32_t u32;

static const u32 MOD = 1000000007u;

static u32 mod_pow(u32 a, u64 e) {
    u64 r = 1, b = a % MOD;
    while (e > 0) {
        if (e & 1) r = r * b % MOD;
        b = b * b % MOD;
        e >>= 1;
    }
    return (u32)r;
}

/* Factorials and inverse factorials */
static u32 *fact, *inv_fact;

static void build_factorials(int N) {
    fact = (u32 *)calloc((size_t)(N + 1), sizeof(u32));
    inv_fact = (u32 *)calloc((size_t)(N + 1), sizeof(u32));
    fact[0] = 1;
    fact[1] = 1;
    for (int i = 2; i <= N; i++)
        fact[i] = (u32)((u64)fact[i - 1] * (u32)i % MOD);
    inv_fact[N] = mod_pow(fact[N], MOD - 2);
    for (int i = N; i >= 1; i--)
        inv_fact[i - 1] = (u32)((u64)inv_fact[i] * (u32)i % MOD);
}

static u32 nCk(int n, int k) {
    if (k < 0 || k > n) return 0;
    return (u32)((u64)fact[n] * inv_fact[k] % MOD * inv_fact[n - k] % MOD);
}

/* Walsh-Hadamard Transform (XOR) in place */
static void fwht_xor(u32 *a, int n) {
    int h = 1;
    while (h < n) {
        int step = h << 1;
        for (int i = 0; i < n; i += step) {
            for (int j = i; j < i + h; j++) {
                u32 x = a[j];
                u32 y = a[j + h];
                u32 u = x + y;
                if (u >= MOD) u -= MOD;
                u32 v;
                if (x >= y) v = x - y;
                else v = x + MOD - y;
                a[j] = u;
                a[j + h] = v;
            }
        }
        h = step;
    }
}

/* Q(X,k) = X*C(X+k+1,k+1) - (k+1)*C(X+k+1,k+2) */
static u32 Q_of(int X, int k) {
    int n1 = X + k + 1;
    u32 c1 = nCk(n1, k + 1);
    u32 c2 = nCk(n1, k + 2);
    i64 val = (i64)X * (i64)c1 - (i64)(k + 1) * (i64)c2;
    val %= (i64)MOD;
    if (val < 0) val += MOD;
    return (u32)val;
}

/* Build single-suit Grundy distribution, padded to length L */
static void build_single_suit(u32 *a, int n, int L) {
    memset(a, 0, (size_t)L * sizeof(u32));
    if (n <= 0) return;
    if (n == 1) {
        a[0] = 2;
        return;
    }

    build_factorials(n);

    u32 pow2_n2 = mod_pow(2, (u64)(n - 2));
    u32 a0 = (u32)((pow2_n2 + 2) % MOD);
    u32 a1 = (u32)((pow2_n2 + (u32)(n - 2)) % MOD);
    a[0] = a0;
    if (n > 1) a[1] = a1;
    if (n > 2) {
        u32 pow2_n3 = mod_pow(2, (u64)(n - 3));
        a[2] = (u32)((pow2_n3 + (u32)(n - 3)) % MOD);
    }

    u32 inv4 = mod_pow(4, MOD - 2);

    /* Odd Grundy values: g = 2k+3, starting at (X0=n-4, k=0) */
    {
        int X0 = n - 4;
        if (3 < n && X0 >= 0) {
            int k = 0;
            int X = X0;
            i64 F = (i64)((mod_pow(2, (u64)(X + 1)) + MOD - 1) % MOD);
            while (1) {
                int g = 2 * k + 3;
                if (g >= n || X < 0) break;
                a[g] = (u32)((F + (i64)Q_of(X, k)) % MOD);

                if (X < 2) break;
                u32 c_xk_1 = nCk(X + k - 1, k);
                u32 c_xk = nCk(X + k, k);
                i64 tmp = F - 2 * (i64)c_xk_1 - (i64)c_xk;
                tmp %= (i64)MOD;
                if (tmp < 0) tmp += MOD;
                tmp = tmp * (i64)inv4 % MOD;
                u32 c_next = nCk(X + k - 1, k + 1);
                F = (2 * tmp - (i64)c_next) % (i64)MOD;
                if (F < 0) F += MOD;
                k++;
                X -= 2;
            }
        }
    }

    /* Even Grundy values: g = 2k+4, starting at (X0=n-5, k=0) */
    {
        int X0 = n - 5;
        if (4 < n && X0 >= 0) {
            int k = 0;
            int X = X0;
            i64 F = (i64)((mod_pow(2, (u64)(X + 1)) + MOD - 1) % MOD);
            while (1) {
                int g = 2 * k + 4;
                if (g >= n || X < 0) break;
                a[g] = (u32)((F + (i64)Q_of(X, k)) % MOD);

                if (X < 2) break;
                u32 c_xk_1 = nCk(X + k - 1, k);
                u32 c_xk = nCk(X + k, k);
                i64 tmp = F - 2 * (i64)c_xk_1 - (i64)c_xk;
                tmp %= (i64)MOD;
                if (tmp < 0) tmp += MOD;
                tmp = tmp * (i64)inv4 % (i64)MOD;
                u32 c_next = nCk(X + k - 1, k + 1);
                F = (2 * tmp - (i64)c_next) % (i64)MOD;
                if (F < 0) F += MOD;
                k++;
                X -= 2;
            }
        }
    }
}

long long p798_native(void) {
    int n = 10000000;
    int s = 10000000;

    if (n == 0) return 1;

    /* L = next power of 2 >= n */
    int L = 1;
    while (L < n) L <<= 1;

    u32 *f = (u32 *)calloc((size_t)L, sizeof(u32));
    build_single_suit(f, n, L);

    /* Free factorials to save memory */
    free(fact);
    free(inv_fact);

    fwht_xor(f, L);

    /* Pointwise exponentiation and sum */
    u64 total = 0;
    for (int i = 0; i < L; i++) {
        total += mod_pow(f[i], (u64)s);
        if ((i & 8191) == 0)
            total %= MOD;
    }
    total %= MOD;

    u32 inv_L = mod_pow((u32)L, MOD - 2);
    u64 result = total * inv_L % MOD;

    free(f);
    return (long long)result;
}
