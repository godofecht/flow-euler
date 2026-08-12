// Project Euler 799: Pentagonal Puzzle
// Find the smallest pentagonal number that can be written as a sum of two
// pentagonal numbers in more than 100 different ways.
//
// Uses: 64-bit Miller-Rabin + Pollard Rho factoring, Cornacchia's algorithm
// for sum of two squares, Gaussian integer representation counting, and a
// segmented sieve to scan m values efficiently.
// Port of the Python reference solver to C.
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef uint64_t u64;
typedef int64_t i64;
typedef __uint128_t u128;

// ---- Modular arithmetic ----

static u64 mulmod(u64 a, u64 b, u64 m) {
    return (u64)((u128)a * b % m);
}

static u64 powmod(u64 base, u64 exp, u64 m) {
    u64 result = 1 % m;
    base %= m;
    while (exp > 0) {
        if (exp & 1) result = mulmod(result, base, m);
        base = mulmod(base, base, m);
        exp >>= 1;
    }
    return result;
}

static u64 gcd_u64(u64 a, u64 b) {
    while (b) { u64 t = a % b; a = b; b = t; }
    return a;
}

// ---- PRNG (deterministic) ----

static u64 rng_state = 88172645463325252ULL;
static u64 xorshift64(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

// ---- Primality (deterministic Miller-Rabin for 64-bit) ----

static int is_prime(u64 n) {
    if (n < 2) return 0;
    static const u64 small[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++) {
        if (n % small[i] == 0) return n == small[i];
    }
    u64 d = n - 1;
    int s = 0;
    while ((d & 1) == 0) { d >>= 1; s++; }
    static const u64 bases[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    for (int i = 0; i < 7; i++) {
        u64 a = bases[i];
        if (a % n == 0) continue;
        u64 x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        int found = 0;
        for (int j = 0; j < s - 1; j++) {
            x = mulmod(x, x, n);
            if (x == n - 1) { found = 1; break; }
        }
        if (!found) return 0;
    }
    return 1;
}

// ---- Pollard Rho ----

static u64 pollard_rho(u64 n) {
    if (n % 2 == 0) return 2;
    if (n % 3 == 0) return 3;
    while (1) {
        u64 c = 1 + xorshift64() % (n - 2);
        u64 x = xorshift64() % n;
        u64 y = x;
        u64 d = 1;
        while (d == 1) {
            x = (mulmod(x, x, n) + c) % n;
            y = (mulmod(y, y, n) + c) % n;
            y = (mulmod(y, y, n) + c) % n;
            d = gcd_u64(x > y ? x - y : y - x, n);
        }
        if (d != n) return d;
    }
}

// ---- Factorize ----

typedef struct { u64 prime; int exp; } PF;

static int factorize(u64 n, PF *out, int max_out) {
    if (n <= 1) return 0;
    u64 stack[128];
    int sp = 0;
    stack[sp++] = n;
    u64 primes[128];
    int np = 0;
    while (sp > 0) {
        u64 x = stack[--sp];
        if (x == 1) continue;
        if (is_prime(x)) {
            primes[np++] = x;
            continue;
        }
        u64 d = pollard_rho(x);
        stack[sp++] = d;
        stack[sp++] = x / d;
    }
    // Sort
    for (int i = 0; i < np; i++)
        for (int j = i + 1; j < np; j++)
            if (primes[j] < primes[i]) {
                u64 tmp = primes[i]; primes[i] = primes[j]; primes[j] = tmp;
            }
    int count = 0;
    for (int i = 0; i < np;) {
        u64 p = primes[i];
        int e = 0;
        while (i < np && primes[i] == p) { e++; i++; }
        if (count < max_out) {
            out[count].prime = p;
            out[count].exp = e;
        }
        count++;
    }
    return count;
}

// ---- Integer sqrt ----

static i64 isqrt_i64(i64 n) {
    if (n <= 0) return 0;
    i64 r = (i64)sqrtl((long double)n);
    while (r > 0 && r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r;
}

// ---- sqrt(-1) mod p (for prime p = 1 mod 4) ----

static u64 sqrt_neg1_mod_p(u64 p) {
    static const u64 bases[] = {2, 3, 5, 6, 7, 10, 11, 13, 17, 19, 23, 29};
    for (int i = 0; i < 12; i++) {
        u64 a = bases[i];
        if (a % p == 0) continue;
        if (powmod(a, (p - 1) / 2, p) == p - 1) {
            return powmod(a, (p - 1) / 4, p);
        }
    }
    u64 a = 2;
    while (1) {
        if (a % p != 0 && powmod(a, (p - 1) / 2, p) == p - 1) {
            return powmod(a, (p - 1) / 4, p);
        }
        a++;
    }
}

// ---- Cornacchia's algorithm: find (a,b) with a^2 + b^2 = p ----

static void cornacchia(u64 p, i64 *out_a, i64 *out_b) {
    u64 t = sqrt_neg1_mod_p(p);
    for (int attempt = 0; attempt < 2; attempt++) {
        if (attempt == 1) t = p - t;
        u64 r0 = p, r1 = t;
        while (r1 * r1 > p) {
            u64 tmp = r0 % r1;
            r0 = r1;
            r1 = tmp;
        }
        i64 a = (i64)r1;
        i64 b2 = (i64)p - a * a;
        i64 b = isqrt_i64(b2);
        if (b * b == b2) {
            *out_a = a;
            *out_b = b;
            return;
        }
    }
    *out_a = 0;
    *out_b = 0;
}

// ---- Gaussian integer multiplication ----

typedef struct { i64 u, v; } GInt;

static GInt gmul(GInt a, GInt b) {
    GInt r;
    r.u = a.u * b.u - a.v * b.v;
    r.v = a.u * b.v + a.v * b.u;
    return r;
}

// ---- Count ways to write P_m as sum of two pentagonal numbers ----

static int cmp_gint(const void *a, const void *b) {
    const GInt *ga = a, *gb = b;
    if (ga->u != gb->u) return ga->u < gb->u ? -1 : 1;
    if (ga->v != gb->v) return ga->v < gb->v ? -1 : 1;
    return 0;
}

static int count_pentagonal_sum_ways(u64 m) {
    u64 x = 6 * m - 1;
    u64 N = x * x + 1;

    PF fac[64];
    int nfac = factorize(N, fac, 64);

    // Build options for each odd prime factor
    GInt options[64][64];
    int opt_sizes[64];
    int n_opts = 0;

    for (int fi = 0; fi < nfac; fi++) {
        u64 p = fac[fi].prime;
        int e = fac[fi].exp;
        if (p == 2) continue;

        i64 a, b;
        cornacchia(p, &a, &b);
        GInt gp = {a, b};
        GInt gc = {a, -b};

        // Compute powers gp^k and gc^k
        GInt pow_gp[64], pow_gc[64];
        pow_gp[0] = (GInt){1, 0};
        pow_gc[0] = (GInt){1, 0};
        for (int k = 1; k <= e; k++) {
            pow_gp[k] = gmul(pow_gp[k - 1], gp);
            pow_gc[k] = gmul(pow_gc[k - 1], gc);
        }

        // opt[k] = gp^k * gc^(e-k) for k=0..e
        for (int k = 0; k <= e; k++) {
            options[n_opts][k] = gmul(pow_gp[k], pow_gc[e - k]);
        }
        opt_sizes[n_opts] = e + 1;
        n_opts++;
    }

    // Cartesian product
    static GInt reps[8192];
    static GInt new_reps[8192];
    int n_reps = 1;
    reps[0] = (GInt){1, 0};

    for (int oi = 0; oi < n_opts; oi++) {
        int n_new = 0;
        for (int ri = 0; ri < n_reps; ri++) {
            for (int ok = 0; ok < opt_sizes[oi]; ok++) {
                if (n_new < 8192) {
                    new_reps[n_new] = gmul(reps[ri], options[oi][ok]);
                    n_new++;
                }
            }
        }
        memcpy(reps, new_reps, n_new * sizeof(GInt));
        n_reps = n_new;
    }

    // Multiply by (1+i) for the factor of 2
    for (int i = 0; i < n_reps; i++) {
        reps[i] = gmul(reps[i], (GInt){1, 1});
    }

    // Filter: both components nonzero, abs, sort u<=v, both u%3==2 and v%3==2
    GInt filtered[8192];
    int n_filtered = 0;
    for (int i = 0; i < n_reps; i++) {
        i64 u = reps[i].u < 0 ? -reps[i].u : reps[i].u;
        i64 v = reps[i].v < 0 ? -reps[i].v : reps[i].v;
        if (u == 0 || v == 0) continue;
        if (u > v) { i64 tmp = u; u = v; v = tmp; }
        if (u % 3 == 2 && v % 3 == 2) {
            filtered[n_filtered].u = u;
            filtered[n_filtered].v = v;
            n_filtered++;
        }
    }

    // Count unique pairs
    qsort(filtered, n_filtered, sizeof(GInt), cmp_gint);
    int seen = 0;
    for (int i = 0; i < n_filtered; i++) {
        if (i == 0 || filtered[i].u != filtered[i-1].u ||
            filtered[i].v != filtered[i-1].v) {
            seen++;
        }
    }
    return seen;
}

// ---- Sieve of Eratosthenes ----

static void sieve_primes(u64 limit, u64 *primes, int *count) {
    char *is_comp = calloc(limit + 1, 1);
    *count = 0;
    for (u64 i = 2; i <= limit; i++) {
        if (!is_comp[i]) {
            primes[(*count)++] = i;
            for (u64 j = i * i; j <= limit; j += i)
                is_comp[j] = 1;
        }
    }
    free(is_comp);
}

// ---- Precompute roots for segmented sieve ----

typedef struct { u64 p, r1, r2; } RootEntry;

static RootEntry g_roots[16384];
static int g_n_roots;
static u64 g_next_prime;

static u64 next_prime_after(u64 n) {
    u64 x = n + 1;
    if (x % 2 == 0) x++;
    while (!is_prime(x)) x += 2;
    return x;
}

static void precompute_roots(u64 prime_limit) {
    u64 primes[20000];
    int n_primes;
    sieve_primes(prime_limit, primes, &n_primes);

    g_n_roots = 0;
    for (int i = 0; i < n_primes; i++) {
        u64 p = primes[i];
        if (p % 4 != 1) continue;
        u64 s = sqrt_neg1_mod_p(p);
        u64 inv6 = powmod(6, p - 2, p);
        u64 r1 = ((1 + s) * inv6) % p;
        u64 r2 = ((1 + p - s) * inv6) % p;
        g_roots[g_n_roots].p = p;
        g_roots[g_n_roots].r1 = r1;
        g_roots[g_n_roots].r2 = r2;
        g_n_roots++;
    }
    g_next_prime = next_prime_after(prime_limit);
}

// ---- Upper bound on multiplicative factor from remaining part ----

static u64 upper_factor_multiplier(u64 rem, u64 min_prime) {
    if (rem <= 1) return 1;
    int k = 0;
    u128 p = min_prime;
    while (p <= rem) {
        p *= min_prime;
        k++;
    }
    return 1ULL << k;
}

// ---- Segmented sieve to find the first m with >100 ways ----

static u64 find_answer(u64 block_size, u64 prime_limit) {
    precompute_roots(prime_limit);

    u64 B = block_size;
    u64 *res = malloc(B * sizeof(u64));
    u64 *prod = malloc(B * sizeof(u64));
    if (!res || !prod) return 0;

    u64 m0 = 1;
    while (1) {
        // Compute residuals: res[i] = ((6*(m0+i)-1)^2 + 1) / 2
        u64 x = 6 * m0 - 1;
        u64 n = x * x + 1;
        for (u64 i = 0; i < B; i++) {
            res[i] = n / 2;
            prod[i] = 1;
            x += 6;
            n = x * x + 1;
        }

        // Divide out small primes using precomputed roots
        for (int ri = 0; ri < g_n_roots; ri++) {
            u64 p = g_roots[ri].p;
            u64 mp = m0 % p;
            for (int r_idx = 0; r_idx < 2; r_idx++) {
                u64 r = (r_idx == 0) ? g_roots[ri].r1 : g_roots[ri].r2;
                u64 idx = (r + p - mp) % p;
                while (idx < B) {
                    u64 val = res[idx];
                    if (val % p == 0) {
                        int e = 0;
                        while (val % p == 0) {
                            val /= p;
                            e++;
                        }
                        res[idx] = val;
                        prod[idx] *= (e + 1);
                    }
                    idx += p;
                }
            }
        }

        // Check candidates
        for (u64 i = 0; i < B; i++) {
            u64 pr = prod[i];
            u64 rem = res[i];
            u64 ub = pr * upper_factor_multiplier(rem, g_next_prime);
            if (ub < 202) continue;

            u64 m = m0 + i;
            int ways = count_pentagonal_sum_ways(m);
            if (ways > 100) {
                free(res);
                free(prod);
                return m * (3 * m - 1) / 2;
            }
        }

        m0 += B;
    }
}

// ---- Entry point ----

long long p799_native(void) {
    return (long long)find_answer(50000, 200000);
}
