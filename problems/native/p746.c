#include <stdint.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
enum { MOD = 1000000007LL };

static i64 modpow(i64 b, i64 e) {
    i64 r = 1; b %= MOD;
    if (b < 0) b += MOD;
    while (e > 0) {
        if (e & 1) r = r * b % MOD;
        b = b * b % MOD;
        e >>= 1;
    }
    return r;
}

static i64 *fact_arr;
static i64 *invfact_arr;
static i64 *inv_arr;
static i64 *pow4_arr;

static void build_factorials(i64 n_max) {
    fact_arr = malloc((size_t)(n_max + 1) * sizeof(i64));
    invfact_arr = malloc((size_t)(n_max + 1) * sizeof(i64));
    fact_arr[0] = 1;
    for (i64 i = 1; i <= n_max; i++) {
        fact_arr[i] = fact_arr[i - 1] * i % MOD;
    }
    invfact_arr[n_max] = modpow(fact_arr[n_max], MOD - 2);
    for (i64 i = n_max; i >= 1; i--) {
        invfact_arr[i - 1] = invfact_arr[i] * i % MOD;
    }
}

static i64 nCk(i64 n, i64 k) {
    if (k < 0 || k > n) return 0;
    return fact_arr[n] * invfact_arr[k] % MOD * invfact_arr[n - k] % MOD;
}

static void precompute_inverses(i64 n_max) {
    inv_arr = malloc((size_t)(n_max + 1) * sizeof(i64));
    memset(inv_arr, 0, (size_t)(n_max + 1) * sizeof(i64));
    inv_arr[1] = 1;
    for (i64 i = 2; i <= n_max; i++) {
        inv_arr[i] = MOD - (MOD / i) * inv_arr[MOD % i] % MOD;
    }
}

static i64 M(i64 n) {
    if (n <= 0) return 0;
    if (n == 1) return 0;
    i64 total = 0;
    i64 nPk = 1;
    for (i64 k = 0; k <= n; k++) {
        if (k > 0) {
            nPk = nPk * (n - (k - 1)) % MOD;
        }
        i64 D;
        if (k == 0) {
            D = 1;
        } else {
            D = (4 * n) % MOD;
            D = D * inv_arr[k] % MOD;
            D = D * nCk(4 * n - 3 * k - 1, k - 1) % MOD;
        }
        i64 rem = 2 * (n - k);
        i64 ways_rest = fact_arr[rem] * fact_arr[rem] % MOD;
        i64 term = nPk;
        term = term * D % MOD;
        term = term * pow4_arr[k] % MOD;
        term = term * ways_rest % MOD;
        if (k & 1) {
            total = (total - term) % MOD;
        } else {
            total = (total + term) % MOD;
        }
    }
    return (2 * total) % MOD;
}

static i64 S(i64 n) {
    i64 acc = 0;
    for (i64 k = 2; k <= n; k++) {
        acc += M(k);
        acc %= MOD;
    }
    return acc;
}

long long p746_native(void) {
    i64 target = 2021;
    i64 max_n = target;
    if (max_n < 10) max_n = 10;
    build_factorials(4 * max_n);
    precompute_inverses(max_n);
    pow4_arr = malloc((size_t)(max_n + 1) * sizeof(i64));
    pow4_arr[0] = 1;
    for (i64 i = 1; i <= max_n; i++) {
        pow4_arr[i] = pow4_arr[i - 1] * 4 % MOD;
    }
    i64 result = S(target) % MOD;
    if (result < 0) result += MOD;
    free(fact_arr);
    free(invfact_arr);
    free(inv_arr);
    free(pow4_arr);
    return result;
}
