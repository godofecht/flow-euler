#include <stdint.h>
#include <string.h>

typedef long long i64;
typedef __int128 i128;

/* ---------- Sieve ---------- */

#define SIEVE_LIMIT 100000
static i64 g_primes[10000];
static int g_prime_count = 0;
static int g_sieve_done = 0;

static void init_sieve(void) {
    if (g_sieve_done) return;
    g_sieve_done = 1;
    static char is_comp[SIEVE_LIMIT + 1];
    memset(is_comp, 0, sizeof(is_comp));
    for (i64 p = 2; p * p <= SIEVE_LIMIT; p++)
        if (!is_comp[p])
            for (i64 i = p * p; i <= SIEVE_LIMIT; i += p)
                is_comp[i] = 1;
    for (i64 i = 2; i <= SIEVE_LIMIT; i++)
        if (!is_comp[i])
            g_primes[g_prime_count++] = i;
}

/* ---------- Factorization ---------- */

typedef struct {
    i64 p[32];
    int e[32];
    int n;
} factors_t;

static factors_t factorize(i64 n) {
    factors_t f;
    f.n = 0;
    if (n <= 1) return f;
    for (int i = 0; i < g_prime_count; i++) {
        i64 p = g_primes[i];
        if (p * p > n) break;
        if (n % p == 0) {
            f.p[f.n] = p;
            f.e[f.n] = 0;
            while (n % p == 0) { n /= p; f.e[f.n]++; }
            f.n++;
        }
    }
    if (n > 1) {
        f.p[f.n] = n;
        f.e[f.n] = 1;
        f.n++;
    }
    return f;
}

static factors_t merge_factors(factors_t a, factors_t b) {
    factors_t r = a;
    for (int i = 0; i < b.n; i++) {
        int found = -1;
        for (int j = 0; j < r.n; j++)
            if (r.p[j] == b.p[i]) { found = j; break; }
        if (found >= 0) r.e[found] += b.e[i];
        else { r.p[r.n] = b.p[i]; r.e[r.n] = b.e[i]; r.n++; }
    }
    return r;
}

static factors_t factorize_t4m1(i64 t) {
    factors_t f1 = factorize(t - 1);
    factors_t f2 = factorize(t + 1);
    factors_t f3 = factorize(t * t + 1);
    factors_t r = merge_factors(f1, f2);
    return merge_factors(r, f3);
}

/* ---------- Number theory ---------- */

static i64 gcd64(i64 a, i64 b) {
    while (b) { i64 t = b; b = a % b; a = t; }
    return a;
}

static i64 lcm64(i64 a, i64 b) {
    return (a / gcd64(a, b)) * b;
}

static i64 powmod128(i64 a, i64 e, i64 mod) {
    i128 r = 1, b = a % mod;
    if (b < 0) b += mod;
    while (e > 0) {
        if (e & 1) r = r * b % (i128)mod;
        b = b * b % (i128)mod;
        e >>= 1;
    }
    return (i64)r;
}

static i64 phi_pp(i64 p, int a) {
    if (a == 0) return 1;
    if (p == 2) return 1LL << (a - 1);
    i64 r = p - 1;
    for (int i = 1; i < a; i++) r *= p;
    return r;
}

static i64 carmichael(i64 p, int exp) {
    if (exp <= 0) return 1;
    if (p == 2) {
        if (exp == 1) return 1;
        if (exp == 2) return 2;
        return 1LL << (exp - 2);
    }
    i64 r = p - 1;
    for (int i = 1; i < exp; i++) r *= p;
    return r;
}

static i64 mult_order(i64 a, i64 p, int exp) {
    i64 mod = 1;
    for (int i = 0; i < exp; i++) mod *= p;
    if (mod == 1) return 1;
    a %= mod;
    i64 lam = carmichael(p, exp);
    i64 q[32]; int qn = 0;
    if (p == 2) {
        if (exp >= 2) q[qn++] = 2;
    } else {
        factors_t pf = factorize(p - 1);
        for (int i = 0; i < pf.n; i++) q[qn++] = pf.p[i];
        if (exp >= 2) q[qn++] = p;
    }
    i64 order = lam;
    for (int i = 0; i < qn; i++) {
        while (order % q[i] == 0 && powmod128(a, order / q[i], mod) == 1)
            order /= q[i];
    }
    return order;
}

/* ---------- Sum phi/ord over divisors ---------- */

#define MAX_PF 20
#define MAX_EX 64

static i64 g_phi[MAX_PF][MAX_EX];
static i64 g_ord[MAX_PF][MAX_EX];
static int g_exps[MAX_PF];

static void dfs(int i, int np, i64 cur_phi, i64 cur_ord, i64 *total) {
    if (i == np) {
        *total += cur_phi / cur_ord;
        return;
    }
    for (int a = 0; a <= g_exps[i]; a++)
        dfs(i + 1, np, cur_phi * g_phi[i][a], lcm64(cur_ord, g_ord[i][a]), total);
}

static i64 sum_phi_over_order(i64 mult, factors_t *fac) {
    int np = fac->n;
    for (int i = 0; i < np; i++) {
        i64 p = fac->p[i];
        int e = fac->e[i];
        g_exps[i] = e;
        for (int a = 0; a <= e; a++) {
            g_phi[i][a] = phi_pp(p, a);
            g_ord[i][a] = (a == 0) ? 1 : mult_order(mult, p, a);
        }
    }
    i64 total = 0;
    dfs(0, np, 1, 1, &total);
    return total;
}

/* ---------- Main ---------- */

long long p913_native(void) {
    init_sieve();
    i128 ans = 0;
    for (int n = 2; n <= 100; n++) {
        for (int m = n; m <= 100; m++) {
            i64 t = (i64)n * m;
            i64 N = t * t * t * t;
            factors_t fac = factorize_t4m1(t);
            i64 mult = (i64)m * m * m * m;
            i64 contrib = sum_phi_over_order(mult, &fac);
            i64 S = (N - 1) - contrib;
            ans += S;
        }
    }
    return (i64)ans;
}
