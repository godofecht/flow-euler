// Project Euler 851: SOP and POS
// R_6(10000!) mod 1e9+7
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

#define MOD 1000000007LL
#define LIMIT 10000

static i64 modpow(i64 base, i64 exp, i64 mod) {
    i64 r = 1; base %= mod; if (base < 0) base += mod;
    while (exp > 0) {
        if (exp & 1) r = (i64)((i128)r * base % mod);
        base = (i64)((i128)base * base % mod);
        exp >>= 1;
    }
    return r;
}

// Sieve
static int is_prime[LIMIT + 1];
static int primes[2000];
static int num_primes;

static void sieve(void) {
    memset(is_prime, 1, sizeof(is_prime));
    is_prime[0] = is_prime[1] = 0;
    for (int p = 2; p * p <= LIMIT; p++) {
        if (is_prime[p]) {
            for (int m = p * p; m <= LIMIT; m += p) is_prime[m] = 0;
        }
    }
    num_primes = 0;
    for (int i = 2; i <= LIMIT; i++) if (is_prime[i]) primes[num_primes++] = i;
}

// Factorial prime exponents: {p: v_p(n!)} for primes p <= n
typedef struct { int p; int e; } PrimeExp;

static PrimeExp prime_exps[2000];
static int num_exps;

static void factorial_prime_exponents(int n) {
    num_exps = 0;
    for (int pi = 0; pi < num_primes; pi++) {
        int p = primes[pi];
        if (p > n) break;
        int e = 0;
        int t = n;
        while (t) { t /= p; e += t; }
        if (e) {
            prime_exps[num_exps].p = p;
            prime_exps[num_exps].e = e;
            num_exps++;
        }
    }
}

// Modular inverses 1..n
static i64 inv_table_arr[LIMIT + 1];

static void compute_inv_table(int n, i64 mod) {
    inv_table_arr[1] = 1;
    for (int i = 2; i <= n; i++)
        inv_table_arr[i] = (mod - (mod / i) * inv_table_arr[mod % i] % mod) % mod;
}

// sigma_s(n) mod mod from prime exponents
static i64 sigma_power_from_exps(i64 s, i64 mod) {
    i64 res = 1;
    for (int i = 0; i < num_exps; i++) {
        int p = prime_exps[i].p;
        int e = prime_exps[i].e;
        i64 ps = modpow(p, s, mod);
        i64 term;
        if (ps == 1) {
            term = (e + 1) % mod;
        } else {
            i64 num = (modpow(ps, e + 1, mod) - 1 + mod) % mod;
            i64 den_inv = modpow((ps - 1 + mod) % mod, mod - 2, mod);
            term = (i64)((i128)num * den_inv % mod);
        }
        res = (i64)((i128)res * term % mod);
    }
    return res;
}

// n mod from prime exponents
static i64 n_mod_from_exps(i64 mod) {
    i64 res = 1;
    for (int i = 0; i < num_exps; i++) {
        res = (i64)((i128)res * modpow(prime_exps[i].p, prime_exps[i].e, mod) % mod);
    }
    return res;
}

// Precompute sigma1[1..n]
static i64 sigma1[LIMIT + 1];

static void precompute_sigma1(int n) {
    memset(sigma1, 0, sizeof(sigma1));
    for (int d = 1; d <= n; d++) {
        for (int m = d; m <= n; m += d) {
            sigma1[m] += d;
            if (sigma1[m] >= MOD) sigma1[m] -= MOD; // Keep bounded (but sigma1 can exceed MOD before mod)
        }
    }
    // Actually sigma1[m] can exceed MOD since we're adding d (up to 10000) many times.
    // Let's just mod after.
    for (int i = 0; i <= n; i++) sigma1[i] %= MOD;
}

// Precompute tau[1..n] mod mod using D(Delta) = E2*Delta
static i64 tau_arr[LIMIT + 1];

static void precompute_tau(int n, i64 mod) {
    precompute_sigma1(n);
    compute_inv_table(n, mod);

    memset(tau_arr, 0, sizeof(tau_arr));
    tau_arr[1] = 1;

    for (int k = 2; k <= n; k++) {
        i64 total = 0;
        for (int m = 1; m < k; m++) {
            total = (total + (i64)((i128)sigma1[m] * tau_arr[k - m] % mod)) % mod;
        }
        i64 val = (i64)((i128)(MOD - 24) % mod * total % mod);
        val = (i64)((i128)val * inv_table_arr[k - 1] % mod);
        tau_arr[k] = val;
    }
}

// 2x2 matrix multiply mod
static void mat_mul(i64 A[4], i64 B[4], i64 mod, i64 R[4]) {
    R[0] = (i64)((i128)A[0] * B[0] % mod + (i128)A[1] * B[2] % mod) % mod;
    R[1] = (i64)((i128)A[0] * B[1] % mod + (i128)A[1] * B[3] % mod) % mod;
    R[2] = (i64)((i128)A[2] * B[0] % mod + (i128)A[3] * B[2] % mod) % mod;
    R[3] = (i64)((i128)A[2] * B[1] % mod + (i128)A[3] * B[3] % mod) % mod;
}

// tau(p^e) mod mod
static i64 tau_prime_power(int p, int e, i64 tau_p, i64 mod) {
    if (e == 0) return 1;
    if (e == 1) return tau_p % mod;

    i64 p11 = modpow(p, 11, mod);
    i64 M[4] = {tau_p % mod, (mod - p11) % mod, 1, 0};
    i64 R[4] = {1, 0, 0, 1};
    int exp = e - 1;
    while (exp) {
        if (exp & 1) {
            i64 tmp[4];
            mat_mul(R, M, mod, tmp);
            memcpy(R, tmp, sizeof(R));
        }
        i64 tmp[4];
        mat_mul(M, M, mod, tmp);
        memcpy(M, tmp, sizeof(M));
        exp >>= 1;
    }
    return (i64)((i128)R[0] * (tau_p % mod) % mod + R[1]) % mod;
}

// tau(n) from prime exponents
static i64 tau_from_exps(i64 mod) {
    i64 res = 1;
    for (int i = 0; i < num_exps; i++) {
        int p = prime_exps[i].p;
        int e = prime_exps[i].e;
        res = (i64)((i128)res * tau_prime_power(p, e, tau_arr[p], mod) % mod);
    }
    return res;
}

// Eisenstein series coefficients
static i64 E_COEFF[13]; // indexed by k=2,4,6,8,10,12

static void init_e_coeff(i64 mod) {
    i64 inv_691 = modpow(691, mod - 2, mod);
    E_COEFF[2] = (mod - 24) % mod;
    E_COEFF[4] = 240 % mod;
    E_COEFF[6] = (mod - 504) % mod;
    E_COEFF[8] = 480 % mod;
    E_COEFF[10] = (mod - 264) % mod;
    E_COEFF[12] = (i64)((i128)(65520 % mod) * inv_691 % mod);
}

// sigma values: sig[1], sig[3], sig[5], sig[7], sig[9], sig[11]
static i64 sig[13]; // indexed by s

static void build_sigma_data(i64 mod) {
    for (int s = 1; s <= 11; s += 2) {
        sig[s] = sigma_power_from_exps(s, mod);
    }
}

// Coefficient of q^n in E_k for n>0
static i64 coeff_Ek(int k, i64 mod) {
    if (k == 2) return (i64)((i128)E_COEFF[2] * sig[1] % mod);
    return (i64)((i128)E_COEFF[k] * sig[k - 1] % mod);
}

// Coefficient of q^n in D^r(E_k), n>0
static i64 coeff_D_Ek(int k, int r, i64 n_pows[], i64 mod) {
    return (i64)((i128)n_pows[r] * coeff_Ek(k, mod) % mod);
}

// Coefficient of q^n in E2^k for k=1..6
static i64 INV2, INV5, INV7, INV24185;

static i64 coeff_E2_pow(int k, i64 n_pows[], i64 tau_n, i64 mod) {
    if (k == 1) return coeff_Ek(2, mod);

    if (k == 2) {
        // E2^2 = E4 + 12 D(E2)
        i64 t1 = coeff_Ek(4, mod);
        i64 t2 = (i64)((i128)12 * coeff_D_Ek(2, 1, n_pows, mod) % mod);
        return (t1 + t2) % mod;
    }

    if (k == 3) {
        // E2^3 = E6 + 9 D(E4) + 72 D^2(E2)
        i64 t1 = coeff_Ek(6, mod);
        i64 t2 = (i64)((i128)9 * coeff_D_Ek(4, 1, n_pows, mod) % mod);
        i64 t3 = (i64)((i128)72 * coeff_D_Ek(2, 2, n_pows, mod) % mod);
        return (t1 + t2 + t3) % mod;
    }

    if (k == 4) {
        // E2^4 = E8 + 8 D(E6) + (216/5) D^2(E4) + 288 D^3(E2)
        i64 c216_5 = (i64)((i128)216 * INV5 % mod);
        i64 t1 = coeff_Ek(8, mod);
        i64 t2 = (i64)((i128)8 * coeff_D_Ek(6, 1, n_pows, mod) % mod);
        i64 t3 = (i64)((i128)c216_5 * coeff_D_Ek(4, 2, n_pows, mod) % mod);
        i64 t4 = (i64)((i128)288 * coeff_D_Ek(2, 3, n_pows, mod) % mod);
        return (t1 + t2 + t3 + t4) % mod;
    }

    if (k == 5) {
        // E2^5 = E10 + (15/2) D(E8) + (240/7) D^2(E6) + 144 D^3(E4) + 864 D^4(E2)
        i64 c15_2 = (i64)((i128)15 * INV2 % mod);
        i64 c240_7 = (i64)((i128)240 * INV7 % mod);
        i64 t1 = coeff_Ek(10, mod);
        i64 t2 = (i64)((i128)c15_2 * coeff_D_Ek(8, 1, n_pows, mod) % mod);
        i64 t3 = (i64)((i128)c240_7 * coeff_D_Ek(6, 2, n_pows, mod) % mod);
        i64 t4 = (i64)((i128)144 * coeff_D_Ek(4, 3, n_pows, mod) % mod);
        i64 t5 = (i64)((i128)864 * coeff_D_Ek(2, 4, n_pows, mod) % mod);
        return (t1 + t2 + t3 + t4 + t5) % mod;
    }

    if (k == 6) {
        // E2^6 = E12 - (4608/24185) Delta + (36/5) D(E10) + 30 D^2(E8)
        //        + (720/7) D^3(E6) + (2592/7) D^4(E4) + (10368/5) D^5(E2)
        i64 c36_5 = (i64)((i128)36 * INV5 % mod);
        i64 c720_7 = (i64)((i128)720 * INV7 % mod);
        i64 c2592_7 = (i64)((i128)2592 * INV7 % mod);
        i64 c10368_5 = (i64)((i128)10368 * INV5 % mod);
        i64 cDelta = (i64)((i128)((mod - 4608) % mod) * INV24185 % mod);
        i64 t1 = coeff_Ek(12, mod);
        i64 t2 = (i64)((i128)cDelta * (tau_n % mod) % mod);
        i64 t3 = (i64)((i128)c36_5 * coeff_D_Ek(10, 1, n_pows, mod) % mod);
        i64 t4 = (i64)((i128)30 * coeff_D_Ek(8, 2, n_pows, mod) % mod);
        i64 t5 = (i64)((i128)c720_7 * coeff_D_Ek(6, 3, n_pows, mod) % mod);
        i64 t6 = (i64)((i128)c2592_7 * coeff_D_Ek(4, 4, n_pows, mod) % mod);
        i64 t7 = (i64)((i128)c10368_5 * coeff_D_Ek(2, 5, n_pows, mod) % mod);
        return (t1 + t2 + t3 + t4 + t5 + t6 + t7) % mod;
    }

    return 0; // should not reach
}

// C(n, k) for small n
static i64 comb_small(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    i64 num = 1, den = 1;
    for (int i = 1; i <= k; i++) {
        num *= (n - (k - i));
        den *= i;
    }
    return num / den;
}

// R_dim(M) mod MOD
static i64 R_dim_at_n(int dim, i64 n_pows[], i64 tau_n, i64 mod) {
    i64 inv12 = modpow(12, mod - 2, mod);
    i64 scale = modpow(inv12, dim, mod);

    i64 s = 0;
    for (int k = 1; k <= dim; k++) {
        i64 ck = comb_small(dim, k);
        i64 term = (i64)((i128)ck * coeff_E2_pow(k, n_pows, tau_n, mod) % mod);
        if (k % 2 == 1) {
            s = (s - term % mod + mod) % mod;
        } else {
            s = (s + term) % mod;
        }
    }
    return (i64)((i128)s * scale % mod);
}

long long p851_native(void) {
    i64 mod = MOD;
    sieve();

    // Initialize constants
    INV2 = (mod + 1) / 2;
    INV5 = modpow(5, mod - 2, mod);
    INV7 = modpow(7, mod - 2, mod);
    INV24185 = modpow(24185, mod - 2, mod);
    init_e_coeff(mod);

    // Precompute tau(p) for primes p <= 10000
    precompute_tau(LIMIT, mod);

    // Main computation: R_6(10000!) mod MOD
    factorial_prime_exponents(LIMIT);
    build_sigma_data(mod);

    i64 nmod = n_mod_from_exps(mod);
    i64 n_pows[6];
    n_pows[0] = 1;
    for (int i = 1; i <= 5; i++)
        n_pows[i] = (i64)((i128)n_pows[i - 1] * nmod % mod);

    i64 tau_fact = tau_from_exps(mod);

    return R_dim_at_n(6, n_pows, tau_fact, mod);
}
