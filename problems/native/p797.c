#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;
typedef uint64_t u64;
typedef uint32_t u32;

static const u32 MOD = 1000000007u;
#define N 10000000

static u32 mod_pow(u32 a, u64 e) {
    u64 r = 1, b = a % MOD;
    while (e > 0) {
        if (e & 1) r = r * b % MOD;
        b = b * b % MOD;
        e >>= 1;
    }
    return (u32)r;
}

/* spf: smallest prime factor, mu: Mobius */
static int *spf;
static signed char *mu;
static u32 *primes;
static int prime_count;

static void linear_sieve(void) {
    spf = (int *)calloc((size_t)(N + 1), sizeof(int));
    mu = (signed char *)calloc((size_t)(N + 1), 1);
    primes = (u32 *)calloc(N / 10 + 1000, sizeof(u32)); /* over-allocate */
    spf[1] = 1;
    mu[1] = 1;
    for (int i = 2; i <= N; i++) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes[prime_count++] = (u32)i;
            mu[i] = -1;
        }
        for (int j = 0; j < prime_count; j++) {
            u32 p = primes[j];
            long ip = (long)i * (long)p;
            if (ip > N) break;
            spf[ip] = (int)p;
            if (i % (int)p == 0) {
                mu[ip] = 0;
                break;
            }
            mu[ip] = -mu[i];
        }
    }
}

/* b[k] = 2^k - 1 mod MOD, invb[k] = inverse of b[k] */
static u32 *b_arr;
static u32 *invb;

static void build_2pow_minus1_and_inv(void) {
    b_arr = (u32 *)calloc((size_t)(N + 1), sizeof(u32));
    invb = (u32 *)calloc((size_t)(N + 1), sizeof(u32));

    u64 pow2 = 1;
    for (int k = 1; k <= N; k++) {
        pow2 = pow2 * 2 % MOD;
        b_arr[k] = (u32)(pow2 - 1); /* pow2 >= 1 always, so pow2-1 >= 0 */
    }

    /* Batch inversion using prefix products */
    invb[0] = 1;
    u64 acc = 1;
    for (int k = 1; k <= N; k++) {
        acc = acc * b_arr[k] % MOD;
        invb[k] = (u32)acc;
    }

    u64 inv_total = mod_pow((u32)acc, MOD - 2);
    for (int k = N; k >= 1; k--) {
        u32 prev_prefix = invb[k - 1];
        invb[k] = (u32)(inv_total * prev_prefix % MOD);
        inv_total = inv_total * b_arr[k] % MOD;
    }
    invb[0] = 0;
}

/* Compute Phi_n(2) mod MOD using the product formula */
static u32 cyclotomic_value_at_2(int n) {
    if (n == 1) return 1;

    /* Extract distinct primes of n */
    int m = n;
    u32 pp[16]; /* at most ~8 distinct primes for n <= 10^7 */
    int npc = 0;
    while (m > 1) {
        int p = spf[m];
        pp[npc++] = (u32)p;
        while (m % p == 0) m /= p;
    }

    /* Enumerate squarefree divisors and their parity */
    /* Use iterative subset enumeration */
    u64 prod[256]; /* 2^8 = 256 max */
    int par[256];
    int cnt = 1;
    prod[0] = 1;
    par[0] = 0;
    for (int i = 0; i < npc; i++) {
        int L = cnt;
        for (int j = 0; j < L; j++) {
            prod[cnt] = prod[j] * pp[i];
            par[cnt] = par[j] ^ 1;
            cnt++;
        }
    }

    u64 res = 1;
    for (int i = 0; i < cnt; i++) {
        int idx = n / (int)prod[i];
        if (par[i] == 0) {
            res = res * b_arr[idx] % MOD;
        } else {
            res = res * invb[idx] % MOD;
        }
    }
    return (u32)res;
}

/* Build T[k] = product_{d|k} (1 + Phi_d(2)), then convert to prefix sums */
static u32 *T;

static void build_T_prefix(void) {
    T = (u32 *)calloc((size_t)(N + 1), sizeof(u32));
    for (int i = 0; i <= N; i++) T[i] = 1;
    T[0] = 0;

    for (int d = 1; d <= N; d++) {
        u32 phi_d = cyclotomic_value_at_2(d);
        u32 fd = phi_d + 1;
        if (fd >= MOD) fd -= MOD;
        for (int m = d; m <= N; m += d) {
            T[m] = (u32)((u64)T[m] * fd % MOD);
        }
    }

    /* Prefix sum in-place */
    u64 run = 0;
    for (int i = 1; i <= N; i++) {
        run += T[i];
        if (run >= MOD) run -= MOD;
        T[i] = (u32)run;
    }
}

long long p797_native(void) {
    linear_sieve();
    build_2pow_minus1_and_inv();
    build_T_prefix();

    /* Free b_arr, invb, spf to save memory */
    free(b_arr);
    free(invb);
    free(spf);

    /* Prefix sums of mu (Mertens function) */
    int *prefix_mu = (int *)calloc((size_t)(N + 1), sizeof(int));
    int run = 0;
    for (int i = 1; i <= N; i++) {
        run += mu[i];
        prefix_mu[i] = run;
    }

    /* Q_n = sum_{d<=n} mu(d) * T_prefix[n/d], grouped by floor division */
    i64 ans = 0;
    int l = 1;
    while (l <= N) {
        int t = N / l;
        int r = N / t;
        int sum_mu = prefix_mu[r] - prefix_mu[l - 1];
        i64 term = (i64)(sum_mu % (int)MOD) * (i64)T[t];
        ans = (ans + term) % (i64)MOD;
        l = r + 1;
    }

    ans %= (i64)MOD;
    if (ans < 0) ans += MOD;

    free(prefix_mu);
    free(mu);
    free(primes);
    free(T);

    return ans;
}
