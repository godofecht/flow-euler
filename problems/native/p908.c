#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Project Euler 908: Clock Sequence II
   Compute C(10^4) mod 1111211113.
   Port of the Python reference solver. */

typedef int64_t i64;
typedef __int128 i128;

#define MOD 1111211113LL
#define N 10000

/* ---------- Sieve primes ---------- */
static int *sieve_primes(int limit, int *count) {
    if (limit < 2) { *count = 0; return NULL; }
    char *is_comp = calloc(limit + 1, 1);
    int *primes = malloc((limit + 1) * sizeof(int));
    int pc = 0;
    for (int p = 2; p <= limit; p++) {
        if (!is_comp[p]) {
            primes[pc++] = p;
            for (long long m = (long long)p * p; m <= limit; m += p)
                is_comp[m] = 1;
        }
    }
    free(is_comp);
    *count = pc;
    return primes;
}

/* ---------- Linear sieve for Mobius ---------- */
static void mobius_upto(int n, int8_t *mu) {
    /* mu[0..n] */
    memset(mu, 0, n + 1);
    int *primes = malloc((n + 1) * sizeof(int));
    int pc = 0;
    char *is_comp = calloc(n + 1, 1);
    mu[1] = 1;
    for (int i = 2; i <= n; i++) {
        if (!is_comp[i]) {
            primes[pc++] = i;
            mu[i] = -1;
        }
        for (int j = 0; j < pc; j++) {
            int p = primes[j];
            long long v = (long long)i * p;
            if (v > n) break;
            is_comp[v] = 1;
            if (i % p == 0) {
                mu[v] = 0;
                break;
            }
            mu[v] = -mu[i];
        }
    }
    free(primes);
    free(is_comp);
}

/* ---------- k(p^e) ---------- */
static i64 k_prime_power(int p, int e) {
    if (e <= 0) return 1;
    if (p == 2) return 1LL << e;
    i64 k = (p + 1) / 2; /* e = 1 */
    for (int exp = 2; exp <= e; exp++) {
        if (exp % 2 == 0)
            k = p * k - (p - 1);
        else
            k = p * k - (p - 1) / 2;
    }
    return k;
}

/* ---------- Generate moduli (m, k(m)) with k(m) <= max_k ---------- */
/* We store pairs as (m, k) where m can be large. Use i64 for m and k. */
typedef struct { i64 m, k; } Pair;

static Pair *pairs_buf = NULL;
static int pairs_count = 0;
static int pairs_cap = 0;

static void add_pair(i64 m, i64 k) {
    if (pairs_count >= pairs_cap) {
        pairs_cap = pairs_cap ? pairs_cap * 2 : 1024;
        pairs_buf = realloc(pairs_buf, pairs_cap * sizeof(Pair));
    }
    pairs_buf[pairs_count].m = m;
    pairs_buf[pairs_count].k = k;
    pairs_count++;
}

/* Options per prime: list of (p^e, k(p^e)) */
typedef struct { i64 m, k; } Opt;
static Opt *options[20001]; /* indexed by prime index */
static int opt_counts[20001];

static void dfs(int start_idx, int num_primes, int *primes,
                i64 m_cur, i64 k_cur, int max_k) {
    add_pair(m_cur, k_cur);

    for (int j = start_idx; j < num_primes; j++) {
        int p = primes[j];
        int oc = opt_counts[j];
        if (oc == 0) continue;
        /* smallest k-factor for this prime */
        if (k_cur * options[j][0].k > max_k) break;

        for (int oi = 0; oi < oc; oi++) {
            i64 mp = options[j][oi].m;
            i64 kp = options[j][oi].k;
            i64 k_new = k_cur * kp;
            if (k_new > max_k) break;
            /* m_cur * mp could overflow i64? m can be up to product of many primes.
               But k <= 10000, and for odd primes k(p)=(p+1)/2, so p <= 2*k-1 ~ 20000.
               The product m is bounded because k is multiplicative and k <= 10000.
               m can be at most around 10^10 or so. Use i128 check. */
            i128 m_new = (i128)m_cur * mp;
            /* m_new fits in i64 since m <= ~10^12 at most for these small k values */
            dfs(j + 1, num_primes, primes, (i64)m_new, k_new, max_k);
        }
    }
}

static void generate_moduli(int max_k, Pair **out_pairs, int *out_count) {
    int num_primes;
    int *primes = sieve_primes(2 * max_k, &num_primes);

    /* Build options for each prime */
    for (int j = 0; j < num_primes; j++) {
        int p = primes[j];
        Opt opts[64];
        int oc = 0;
        if (p == 2) {
            i64 m = 2, k = 2;
            while (k <= max_k) {
                opts[oc].m = m; opts[oc].k = k; oc++;
                m <<= 1;
                k <<= 1;
            }
        } else {
            i64 m = p, k = (p + 1) / 2;
            int e = 1;
            while (k <= max_k) {
                opts[oc].m = m; opts[oc].k = k; oc++;
                e++;
                m *= p;
                if (e % 2 == 0)
                    k = p * k - (p - 1);
                else
                    k = p * k - (p - 1) / 2;
            }
        }
        options[j] = malloc(oc * sizeof(Opt));
        memcpy(options[j], opts, oc * sizeof(Opt));
        opt_counts[j] = oc;
    }

    pairs_buf = NULL;
    pairs_count = 0;
    pairs_cap = 0;

    dfs(0, num_primes, primes, 1, 1, max_k);

    /* Cleanup options */
    for (int j = 0; j < num_primes; j++) free(options[j]);
    free(primes);

    *out_pairs = pairs_buf;
    *out_count = pairs_count;
}

/* ---------- Modular inverses ---------- */
static i64 *prepare_inverses(int n, i64 mod) {
    i64 *inv = calloc(n + 1, sizeof(i64));
    inv[1] = 1;
    for (int i = 2; i <= n; i++)
        inv[i] = (mod - (mod / i) * inv[mod % i] % mod) % mod;
    return inv;
}

/* ---------- Compute B ---------- */
static void compute_B(int max_period, i64 mod, i64 *B) {
    Pair *moduli;
    int num_moduli;
    generate_moduli(max_period, &moduli, &num_moduli);

    i64 *inv = prepare_inverses(max_period, mod);

    memset(B, 0, (max_period + 1) * sizeof(i64));

    for (int idx = 0; idx < num_moduli; idx++) {
        i64 m = moduli[idx].m;
        i64 k = moduli[idx].k;
        if (k > max_period) continue;
        i64 n = m - k;
        if (n < 0) continue;
        i64 rmax = max_period - k;
        if (rmax < 0) continue;
        if (n < rmax) rmax = n;

        /* r = 0 */
        i64 idx0 = k;
        B[idx0] += 1;
        if (B[idx0] >= mod) B[idx0] -= mod;

        i64 c = 1;
        for (i64 r = 1; r <= rmax; r++) {
            /* c = c * (n - r + 1) * inv[r] % mod */
            c = (i64)((i128)c * (n - r + 1) % mod);
            c = (i64)((i128)c * inv[r] % mod);
            i64 bidx = idx0 + r;
            B[bidx] += c;
            if (B[bidx] >= mod) B[bidx] -= mod;
        }
    }

    free(inv);
    free(moduli);
}

/* ---------- Compute A from B via Mobius ---------- */
static void compute_A_from_B(i64 *B, int8_t *mu, i64 mod, i64 *A, int n) {
    memset(A, 0, (n + 1) * sizeof(i64));
    for (int d = 1; d <= n; d++) {
        int8_t md = mu[d];
        if (md == 0) continue;
        if (md == 1) {
            for (int q = 1; q <= n / d; q++) {
                int p = d * q;
                A[p] += B[q];
                if (A[p] >= mod) A[p] -= mod;
            }
        } else { /* md == -1 */
            for (int q = 1; q <= n / d; q++) {
                int p = d * q;
                A[p] -= B[q];
                if (A[p] < 0) A[p] += mod;
            }
        }
    }
}

long long p908_native(void) {
    i64 *B = calloc(N + 1, sizeof(i64));
    compute_B(N, MOD, B);

    int8_t *mu = calloc(N + 1, 1);
    mobius_upto(N, mu);

    i64 *A = calloc(N + 1, sizeof(i64));
    compute_A_from_B(B, mu, MOD, A, N);

    /* Prefix sum */
    i64 s = 0;
    for (int i = 1; i <= N; i++) {
        s += A[i];
        s %= MOD;
    }

    i64 result = s % MOD;

    free(B); free(mu); free(A);
    return result;
}
