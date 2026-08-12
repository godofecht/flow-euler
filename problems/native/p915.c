// Project Euler 915: Giant GCDs
// s(1)=1, s(n+1) = (s(n)-1)^3 + 2
// T(N) = sum_{a=1..N} sum_{b=1..N} gcd(s(s(a)), s(s(b)))
// Compute T(10^8) mod 123456789.
// Uses cycle detection for s(n) mod M, summatory totient via sieve + memoized recursion.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned long long u64;
typedef __int128 i128;

#define MOD 123456789ULL

static u64 mulmod(u64 a, u64 b, u64 m) {
    return (u64)((i128)a * b % m);
}

static u64 f_mod(u64 x, u64 m) {
    /* (x - 1) % m with proper handling of unsigned underflow.
       In Python, (-1) % m = m-1. In C, (u64)(0-1) = 2^64-1, wrong. */
    u64 y = (x + m - 1) % m;
    u64 y2 = mulmod(y, y, m);
    u64 y3 = mulmod(y2, y, m);
    return (y3 + 2) % m;
}

/* Floyd cycle detection for s(n) mod m, starting from x0=0 (s(0)=0).
   Returns (mu, lam) where the sequence becomes periodic. */
static void cycle_info(u64 m, u64 *mu_out, u64 *lam_out) {
    u64 tortoise = f_mod(0, m);
    u64 hare = f_mod(f_mod(0, m), m);
    while (tortoise != hare) {
        tortoise = f_mod(tortoise, m);
        hare = f_mod(f_mod(hare, m), m);
    }
    u64 mu = 0;
    tortoise = 0;
    while (tortoise != hare) {
        tortoise = f_mod(tortoise, m);
        hare = f_mod(hare, m);
        mu++;
    }
    u64 lam = 1;
    hare = f_mod(tortoise, m);
    while (tortoise != hare) {
        hare = f_mod(hare, m);
        lam++;
    }
    *mu_out = mu;
    *lam_out = lam;
}

static u64 *build_s_mod(u64 m, u64 length) {
    u64 *arr = malloc((length + 1) * sizeof(u64));
    u64 x = 0;
    for (u64 n = 1; n <= length; n++) {
        x = f_mod(x, m);
        arr[n] = x;
    }
    return arr;
}

/* Sieve phi prefix up to n */
static u64 *sieve_phi_prefix(u64 n) {
    u64 *phi = malloc((n + 1) * sizeof(u64));
    for (u64 i = 0; i <= n; i++) phi[i] = i;
    for (u64 i = 2; i <= n; i++) {
        if (phi[i] == i) {
            for (u64 j = i; j <= n; j += i)
                phi[j] -= phi[j] / i;
        }
    }
    u64 *pref = malloc((n + 1) * sizeof(u64));
    u64 s = 0;
    for (u64 i = 1; i <= n; i++) {
        s += phi[i];
        pref[i] = s;
    }
    free(phi);
    return pref;
}

/* ---- Solver ---- */

static u64 *s_mod_MOD_arr;
static u64 muM, lamM;
static u64 *s_mod_lam_arr;
static u64 muP, lamP;
static u64 *small_exact_arr;
static u64 n_small_max;
static u64 muM_mod;

static u64 s_index_modMOD(u64 k) {
    if (k <= muM + lamM)
        return s_mod_MOD_arr[k];
    u64 k2 = muM + ((k - muM) % lamM);
    return s_mod_MOD_arr[k2];
}

static u64 s_n_mod_lamM(u64 n_) {
    if (n_ <= muP + lamP)
        return s_mod_lam_arr[n_];
    u64 n2 = muP + ((n_ - muP) % lamP);
    return s_mod_lam_arr[n2];
}

static u64 s2_mod(u64 n_) {
    if (n_ <= n_small_max)
        return s_index_modMOD(small_exact_arr[n_]);
    u64 k_mod = s_n_mod_lamM(n_);
    /* (k_mod - muM_mod) % lamM with proper handling of unsigned underflow */
    u64 diff;
    if (k_mod >= muM_mod) diff = (k_mod - muM_mod) % lamM;
    else diff = lamM - ((muM_mod - k_mod) % lamM);
    if (diff == lamM) diff = 0;
    u64 idx = muM + diff;
    return s_mod_MOD_arr[idx];
}

/* ---- Summatory totient with memoization ---- */

#define BASE 2000000ULL
static u64 *phi_prefix_arr;
static u64 N_global;

typedef struct { u64 key; u64 val; } Entry;

#define MEMO_CAP (1 << 20)
static Entry *memo_table;

static u64 mix64(u64 x) {
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    x ^= x >> 31; return x;
}

static u64 phi_sum_memo(u64 n_);

static u64 phi_sum_memo(u64 n_) {
    if (n_ <= BASE)
        return phi_prefix_arr[n_];

    u64 h = mix64(n_) & (MEMO_CAP - 1);
    while (memo_table[h].key != 0) {
        if (memo_table[h].key == n_)
            return memo_table[h].val;
        h = (h + 1) & (MEMO_CAP - 1);
    }

    // Compute exact phi_sum (fits in u64 for n <= 10^8)
    // phi_sum(n) = n*(n+1)/2 - sum_{l=2..n} (r-l+1) * phi_sum(n/l)
    // For n=10^8, n*(n+1)/2 ~ 5*10^15, fits in u64
    u64 res = (n_ * (n_ + 1)) / 2;
    u64 l = 2;
    while (l <= n_) {
        u64 q = n_ / l;
        u64 r = n_ / q;
        u64 sub = phi_sum_memo(q);
        u64 cnt = r - l + 1;
        // res -= cnt * sub  (exact, no mod)
        // cnt * sub could overflow u64? cnt <= 10^8, sub <= 5*10^15
        // cnt * sub could be up to 5*10^23, which overflows u64!
        // Use __int128 for the multiplication
        i128 val = (i128)cnt * sub;
        res = (u64)((i128)res - val);
        l = r + 1;
    }
    memo_table[h].key = n_;
    memo_table[h].val = res;
    return res;
}

static u64 coprime_pairs(u64 m_) {
    return (2 * (phi_sum_memo(m_) % MOD) - 1 + MOD) % MOD;
}

long long p915_native(void) {
    N_global = 100000000ULL;

    /* Periodicity of s(n) mod MOD */
    cycle_info(MOD, &muM, &lamM);
    s_mod_MOD_arr = build_s_mod(MOD, muM + lamM);

    /* Periodicity of s(n) mod lamM */
    cycle_info(lamM, &muP, &lamP);
    s_mod_lam_arr = build_s_mod(lamM, muP + lamP);

    /* Compute exact s(n) for small n until s(n) > muM */
    small_exact_arr = malloc(sizeof(u64) * 4);
    u64 scap = 4;
    small_exact_arr[0] = 0;
    u64 x = 0;
    u64 n = 0;
    while (1) {
        n++;
        /* x = (x-1)^3 + 2, exact (no mod) */
        i128 xm1 = (i128)x - 1;
        i128 xc = xm1 * xm1 * xm1 + 2;
        x = (u64)xc;
        if (n >= scap) {
            scap *= 2;
            small_exact_arr = realloc(small_exact_arr, scap * sizeof(u64));
        }
        small_exact_arr[n] = x;
        if (x > muM) break;
    }
    n_small_max = n - 1;
    muM_mod = muM % lamM;

    /* s2_mod becomes periodic once s(n) mod lamM is in its cycle */
    u64 start = muP;
    if (n_small_max + 1 > start) start = n_small_max + 1;
    if (1 > start) start = 1;
    u64 period = lamP;

    u64 *period_vals = malloc(period * sizeof(u64));
    for (u64 i = 0; i < period; i++)
        period_vals[i] = s2_mod(start + i) % MOD;

    /* Prefix sums of s2 for O(1) range sums */
    u64 *small_prefix = calloc(start, sizeof(u64));
    u64 acc = 0;
    for (u64 i = 1; i < start; i++) {
        acc = (acc + s2_mod(i)) % MOD;
        small_prefix[i] = acc;
    }

    u64 *period_prefix = malloc((period + 1) * sizeof(u64));
    u64 accp = 0;
    for (u64 i = 0; i < period; i++) {
        accp = (accp + period_vals[i]) % MOD;
        period_prefix[i + 1] = accp;
    }
    u64 period_sum = period_prefix[period];

    /* prefix_s2(n) = sum_{i=1..n} s(s(i)) mod MOD */
    /* Inlined in the main loop below */

    /* Summatory totient */
    phi_prefix_arr = sieve_phi_prefix(BASE);
    memo_table = calloc(MEMO_CAP, sizeof(Entry));

    /* Block over d where floor(N/d) is constant */
    u64 ans = 0;
    u64 l = 1;
    while (l <= N_global) {
        u64 q = N_global / l;
        u64 r = N_global / q;

        /* sum_s2 = prefix_s2(r) - prefix_s2(l-1) */
        u64 pr, pl;
        /* prefix_s2(r) */
        if (r == 0) pr = 0;
        else if (r < start) pr = small_prefix[r];
        else {
            u64 base = (start > 0) ? small_prefix[start - 1] : 0;
            u64 t = r - (start - 1);
            u64 full = t / period;
            u64 rem = t % period;
            pr = (base + mulmod(full, period_sum, MOD) + period_prefix[rem]) % MOD;
        }
        /* prefix_s2(l-1) */
        u64 lm1 = l - 1;
        if (lm1 == 0) pl = 0;
        else if (lm1 < start) pl = small_prefix[lm1];
        else {
            u64 base = (start > 0) ? small_prefix[start - 1] : 0;
            u64 t = lm1 - (start - 1);
            u64 full = t / period;
            u64 rem = t % period;
            pl = (base + mulmod(full, period_sum, MOD) + period_prefix[rem]) % MOD;
        }
        u64 sum_s2 = (pr + MOD - pl) % MOD;
        ans = (ans + mulmod(sum_s2, coprime_pairs(q), MOD)) % MOD;
        l = r + 1;
    }

    return (long long)(ans % MOD);
}
