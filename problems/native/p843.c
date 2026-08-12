// Project Euler 843: Periodic Circles
// S(100) = sum of all possible eventual periods for ring sizes 3..100
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

// ==================== 256-bit GF(2) polynomials ====================

typedef struct { uint64_t w[4]; } Poly;

static int poly_deg(const Poly *p) {
    for (int i = 3; i >= 0; i--)
        if (p->w[i]) return 63 - __builtin_clzll(p->w[i]) + 64 * i;
    return -1;
}

static int poly_is_zero(const Poly *p) {
    return !(p->w[0] | p->w[1] | p->w[2] | p->w[3]);
}

static int poly_is_one(const Poly *p) {
    return p->w[0] == 1 && !(p->w[1] | p->w[2] | p->w[3]);
}

static void poly_zero(Poly *p) {
    p->w[0] = p->w[1] = p->w[2] = p->w[3] = 0;
}

static void poly_one(Poly *p) {
    p->w[0] = 1; p->w[1] = p->w[2] = p->w[3] = 0;
}

static void poly_copy(Poly *dst, const Poly *src) {
    dst->w[0] = src->w[0]; dst->w[1] = src->w[1];
    dst->w[2] = src->w[2]; dst->w[3] = src->w[3];
}

static int poly_equal(const Poly *a, const Poly *b) {
    return a->w[0] == b->w[0] && a->w[1] == b->w[1] &&
           a->w[2] == b->w[2] && a->w[3] == b->w[3];
}

static int poly_get_bit(const Poly *p, int i) {
    return (int)((p->w[i >> 6] >> (i & 63)) & 1);
}

static void poly_set_bit(Poly *p, int i) {
    p->w[i >> 6] |= (1ULL << (i & 63));
}

static void poly_xor(Poly *dst, const Poly *src) {
    dst->w[0] ^= src->w[0]; dst->w[1] ^= src->w[1];
    dst->w[2] ^= src->w[2]; dst->w[3] ^= src->w[3];
}

// dst = src << shift
static void poly_shl(Poly *dst, const Poly *src, int shift) {
    int ws = shift >> 6, bs = shift & 63;
    poly_zero(dst);
    if (bs == 0) {
        for (int i = 3; i >= 0; i--) {
            int di = i + ws;
            if (di >= 0 && di < 4) dst->w[di] = src->w[i];
        }
    } else {
        for (int i = 3; i >= 0; i--) {
            int di = i + ws;
            if (di >= 0 && di < 4)
                dst->w[di] |= src->w[i] << bs;
            if (di + 1 >= 0 && di + 1 < 4)
                dst->w[di + 1] |= src->w[i] >> (64 - bs);
        }
    }
}

// dst = a * b (carry-less multiply)
static void poly_mul(Poly *dst, const Poly *a, const Poly *b) {
    poly_zero(dst);
    Poly temp;
    int db = poly_deg(b);
    for (int i = 0; i <= db; i++) {
        if (poly_get_bit(b, i)) {
            poly_shl(&temp, a, i);
            poly_xor(dst, &temp);
        }
    }
}

// dst = a mod mod (mod is monic)
static void poly_mod(Poly *dst, const Poly *a, const Poly *mod) {
    Poly r;
    poly_copy(&r, a);
    int md = poly_deg(mod);
    while (!poly_is_zero(&r) && poly_deg(&r) >= md) {
        int shift = poly_deg(&r) - md;
        Poly shifted;
        poly_shl(&shifted, mod, shift);
        poly_xor(&r, &shifted);
    }
    poly_copy(dst, &r);
}

// dst = (a * b) mod mod
static void poly_mul_mod(Poly *dst, const Poly *a, const Poly *b, const Poly *mod) {
    Poly prod;
    poly_mul(&prod, a, b);
    poly_mod(dst, &prod, mod);
}

// dst = a^2 (insert zeros between bits)
static void poly_square(Poly *dst, const Poly *a) {
    poly_zero(dst);
    int da = poly_deg(a);
    for (int i = 0; i <= da; i++)
        if (poly_get_bit(a, i))
            poly_set_bit(dst, 2 * i);
}

// dst = a^2 mod mod
static void poly_square_mod(Poly *dst, const Poly *a, const Poly *mod) {
    Poly sq;
    poly_square(&sq, a);
    poly_mod(dst, &sq, mod);
}

// dst = a^exp mod mod
static void poly_pow_mod(Poly *dst, const Poly *a, i64 exp, const Poly *mod) {
    Poly res, base;
    poly_one(&res);
    poly_mod(&base, a, mod);
    while (exp > 0) {
        if (exp & 1) {
            Poly tmp;
            poly_mul_mod(&tmp, &res, &base, mod);
            poly_copy(&res, &tmp);
        }
        exp >>= 1;
        if (exp) {
            Poly tmp;
            poly_square_mod(&tmp, &base, mod);
            poly_copy(&base, &tmp);
        }
    }
    poly_copy(dst, &res);
}

// gcd(a, b) over GF(2)
static void poly_gcd(Poly *dst, const Poly *a, const Poly *b) {
    Poly x, y, r;
    poly_copy(&x, a);
    poly_copy(&y, b);
    while (!poly_is_zero(&y)) {
        poly_mod(&r, &x, &y);
        poly_copy(&x, &y);
        poly_copy(&y, &r);
    }
    poly_copy(dst, &x);
}

// exact division a / b (b nonzero, monic)
static void poly_div_exact(Poly *dst, const Poly *a, const Poly *b) {
    Poly r, q;
    poly_copy(&r, a);
    poly_zero(&q);
    int db = poly_deg(b);
    while (!poly_is_zero(&r) && poly_deg(&r) >= db) {
        int shift = poly_deg(&r) - db;
        poly_set_bit(&q, shift);
        Poly shifted;
        poly_shl(&shifted, b, shift);
        poly_xor(&r, &shifted);
    }
    poly_copy(dst, &q);
}

// ==================== Small integer factoring ====================

#define SIEVE_LIMIT 2000000
static char sieve_arr[SIEVE_LIMIT + 1];
static int small_primes[200000];
static int num_small_primes;

static void init_sieve(void) {
    memset(sieve_arr, 1, sizeof(sieve_arr));
    sieve_arr[0] = sieve_arr[1] = 0;
    for (int p = 2; (i64)p * p <= SIEVE_LIMIT; p++) {
        if (sieve_arr[p]) {
            for (int m = p * p; m <= SIEVE_LIMIT; m += p)
                sieve_arr[m] = 0;
        }
    }
    num_small_primes = 0;
    for (int i = 2; i <= SIEVE_LIMIT; i++)
        if (sieve_arr[i]) small_primes[num_small_primes++] = i;
}

// Factor n into prime factors with multiplicity. Returns count.
static int factor_small(i64 n, i64 *factors) {
    int cnt = 0;
    for (int i = 0; i < num_small_primes; i++) {
        int p = small_primes[i];
        if ((i64)p * p > n) break;
        while (n % p == 0) {
            factors[cnt++] = p;
            n /= p;
        }
    }
    if (n > 1) factors[cnt++] = n;
    return cnt;
}

// ==================== Berlekamp factorization ====================

// Nullspace basis of (Q - I) where Q is Frobenius a -> a^2 mod f
// Returns basis vectors in basis[], returns count.
static int berlekamp_nullspace(const Poly *f, Poly *basis) {
    int n = poly_deg(f);
    if (n <= 0) return 0;

    // Build Q - I matrix: rows[i] is row i (n bits, use __int128)
    // Column j is x^(2j) mod f
    __int128 rows[128];
    memset(rows, 0, sizeof(rows));

    for (int j = 0; j < n; j++) {
        Poly xpow;
        poly_zero(&xpow);
        poly_set_bit(&xpow, 2 * j);
        Poly col;
        poly_mod(&col, &xpow, f);
        int i = 0;
        while (!poly_is_zero(&col)) {
            if (col.w[0] & 1)
                rows[i] |= ((__int128)1 << j);
            // shift col right by 1
            col.w[0] = (col.w[0] >> 1) | ((col.w[1] & 1) << 63);
            col.w[1] = (col.w[1] >> 1) | ((col.w[2] & 1) << 63);
            col.w[2] = (col.w[2] >> 1) | ((col.w[3] & 1) << 63);
            col.w[3] >>= 1;
            i++;
        }
    }

    // Q - I
    for (int i = 0; i < n; i++)
        rows[i] ^= ((__int128)1 << i);

    // Gauss-Jordan over GF(2)
    int pivot_cols[128];
    int pivot_row_for_col[128];
    int num_pivots = 0;
    int r = 0;
    for (int c = 0; c < n && r < n; c++) {
        int pivot = -1;
        for (int i = r; i < n; i++) {
            if ((rows[i] >> c) & 1) { pivot = i; break; }
        }
        if (pivot < 0) continue;
        __int128 tmp = rows[r]; rows[r] = rows[pivot]; rows[pivot] = tmp;
        __int128 pv = rows[r];
        for (int i = 0; i < n; i++) {
            if (i != r && ((rows[i] >> c) & 1))
                rows[i] ^= pv;
        }
        pivot_cols[num_pivots] = c;
        pivot_row_for_col[c] = r;
        num_pivots++;
        r++;
    }

    // Free columns = non-pivot columns
    int free_cols[128];
    int num_free = 0;
    for (int c = 0; c < n; c++) {
        int is_pivot = 0;
        for (int j = 0; j < num_pivots; j++)
            if (pivot_cols[j] == c) { is_pivot = 1; break; }
        if (!is_pivot) free_cols[num_free++] = c;
    }

    for (int fi = 0; fi < num_free; fi++) {
        int fc = free_cols[fi];
        Poly v;
        poly_zero(&v);
        poly_set_bit(&v, fc);
        for (int j = 0; j < num_pivots; j++) {
            int pc = pivot_cols[j];
            __int128 row = rows[pivot_row_for_col[pc]];
            if ((row >> fc) & 1)
                poly_set_bit(&v, pc);
        }
        poly_copy(&basis[fi], &v);
    }
    return num_free;
}

// Factor a square-free monic polynomial into irreducibles.
// Returns factors in out[], returns count.
static int factor_squarefree(const Poly *f, Poly *out);

static int factor_squarefree(const Poly *f, Poly *out) {
    int n = poly_deg(f);
    if (n <= 0) return 0;
    if (n == 1) { poly_copy(&out[0], f); return 1; }

    Poly basis[128];
    int nbasis = berlekamp_nullspace(f, basis);
    if (nbasis <= 1) { poly_copy(&out[0], f); return 1; }

    // Try all non-empty subsets of basis[1..nbasis-1]
    int total_subsets = 1 << (nbasis - 1);
    for (int mask = 1; mask < total_subsets; mask++) {
        Poly comb;
        poly_zero(&comb);
        for (int j = 0; j < nbasis - 1; j++) {
            if (mask & (1 << j))
                poly_xor(&comb, &basis[j + 1]);
        }
        if (poly_is_zero(&comb) || poly_is_one(&comb))
            continue;

        Poly g;
        poly_gcd(&g, f, &comb);
        int dg = poly_deg(&g);
        if (dg > 0 && dg < n) {
            Poly h;
            poly_div_exact(&h, f, &g);
            int ng = factor_squarefree(&g, out);
            Poly out2[128];
            int nh = factor_squarefree(&h, out2);
            for (int i = 0; i < nh; i++)
                poly_copy(&out[ng + i], &out2[i]);
            return ng + nh;
        }

        // Try comb ^ 1
        Poly comb1;
        poly_copy(&comb1, &comb);
        comb1.w[0] ^= 1;  // comb ^ 1
        poly_gcd(&g, f, &comb1);
        dg = poly_deg(&g);
        if (dg > 0 && dg < n) {
            Poly h;
            poly_div_exact(&h, f, &g);
            int ng = factor_squarefree(&g, out);
            Poly out2[128];
            int nh = factor_squarefree(&h, out2);
            for (int i = 0; i < nh; i++)
                poly_copy(&out[ng + i], &out2[i]);
            return ng + nh;
        }
    }

    // Should not reach here for valid square-free polynomials
    poly_copy(&out[0], f);
    return 1;
}

// ==================== Factor cache for x^m + 1 ====================

static Poly factor_cache[100][128];
static int factor_cache_count[100];
static char factor_cache_done[100];

static int irreducible_factors_xm_plus_1(int m, Poly *out) {
    if (factor_cache_done[m]) {
        for (int i = 0; i < factor_cache_count[m]; i++)
            poly_copy(&out[i], &factor_cache[m][i]);
        return factor_cache_count[m];
    }
    Poly f;
    poly_zero(&f);
    poly_set_bit(&f, m);
    poly_set_bit(&f, 0);  // x^m + 1
    int cnt = factor_squarefree(&f, out);
    for (int i = 0; i < cnt; i++)
        poly_copy(&factor_cache[m][i], &out[i]);
    factor_cache_count[m] = cnt;
    factor_cache_done[m] = 1;
    return cnt;
}

// ==================== Order computations ====================

// Smallest k > 0 such that a^(2^k) == a in GF(2)[x]/(mod_irred)
static int frobenius_orbit_degree(const Poly *a, const Poly *mod_irred) {
    Poly t, sq;
    poly_copy(&t, a);
    int d = poly_deg(mod_irred);
    for (int k = 1; k <= d; k++) {
        poly_square_mod(&sq, &t, mod_irred);
        poly_copy(&t, &sq);
        if (poly_equal(&t, a)) return k;
    }
    return d;
}

// Mersenne factors cache
static i64 mersenne_cache[42][64];
static int mersenne_cache_count[42];
static char mersenne_cache_done[42];

static int mersenne_factors(int k, i64 *out) {
    if (mersenne_cache_done[k]) {
        for (int i = 0; i < mersenne_cache_count[k]; i++)
            out[i] = mersenne_cache[k][i];
        return mersenne_cache_count[k];
    }
    i64 m = (1LL << k) - 1;
    int cnt = factor_small(m, out);
    // Sort
    for (int i = 0; i < cnt - 1; i++)
        for (int j = i + 1; j < cnt; j++)
            if (out[j] < out[i]) { i64 tmp = out[i]; out[i] = out[j]; out[j] = tmp; }
    for (int i = 0; i < cnt; i++)
        mersenne_cache[k][i] = out[i];
    mersenne_cache_count[k] = cnt;
    mersenne_cache_done[k] = 1;
    return cnt;
}

// Order of nonzero a modulo irreducible mod_irred
static i64 multiplicative_order(const Poly *a, const Poly *mod_irred) {
    if (poly_is_zero(a)) return 0;

    int k = frobenius_orbit_degree(a, mod_irred);
    i64 group_order = (1LL << k) - 1;
    if (group_order == 1) return 1;

    i64 order = group_order;
    i64 facs[64];
    int nfac = mersenne_factors(k, facs);
    // Get unique prime factors
    i64 unique_primes[64];
    int nunique = 0;
    for (int i = 0; i < nfac; i++) {
        if (i == 0 || facs[i] != facs[i - 1])
            unique_primes[nunique++] = facs[i];
    }
    for (int i = 0; i < nunique; i++) {
        i64 p = unique_primes[i];
        while (order % p == 0) {
            i64 cand = order / p;
            Poly result;
            poly_pow_mod(&result, a, cand, mod_irred);
            if (poly_is_one(&result))
                order = cand;
            else
                break;
        }
    }
    return order;
}

// For modulus p^t, t=1..max_exp, find max v2(order_t / base_odd_order)
static int max_two_lift_exponent(const Poly *base_poly, i64 base_odd_order,
                                  const Poly *p, int max_exp) {
    if (max_exp <= 1) return 0;

    i64 order = base_odd_order;
    Poly mod, tmp;
    poly_copy(&mod, p);
    for (int t = 2; t <= max_exp; t++) {
        poly_mul(&tmp, &mod, p);
        poly_copy(&mod, &tmp);
        Poly result;
        poly_pow_mod(&result, base_poly, order, &mod);
        while (!poly_is_one(&result)) {
            order *= 2;
            poly_pow_mod(&result, base_poly, order, &mod);
        }
    }

    i64 ratio = order / base_odd_order;
    // ratio is a power of two; return log2(ratio)
    int exp = 0;
    while (ratio > 1) { ratio >>= 1; exp++; }
    return exp;
}

// ==================== Integer GCD/LCM ====================

static i64 igcd(i64 a, i64 b) {
    while (b) { i64 t = a % b; a = b; b = t; }
    return a;
}

static i64 ilcm(i64 a, i64 b) {
    if (a == 0 || b == 0) return 0;
    i128 r = (i128)(a / igcd(a, b)) * b;
    return (i64)r;  // Should fit in i64 for our problem sizes
}

// ==================== Period enumeration ====================

// DP entry: (odd_lcm, max_two_exp)
typedef struct { i64 odd_lcm; int smax; } DPEntry;

static i64 periods_set[100000];
static int periods_set_count;

// Insert into periods set (sorted, no duplicates)
static void periods_insert(i64 val) {
    // Linear search (set is small)
    for (int i = 0; i < periods_set_count; i++)
        if (periods_set[i] == val) return;
    periods_set[periods_set_count++] = val;
}

static void periods_for_n(int n, i64 *out_periods, int *out_count) {
    // n = 2^a * m with m odd
    int m = n, a = 0;
    while (m % 2 == 0) { m /= 2; a++; }
    int max_exp = 1 << a;

    // g(x) = x + x^(n-1)
    Poly g_poly;
    poly_zero(&g_poly);
    poly_set_bit(&g_poly, 1);
    poly_set_bit(&g_poly, n - 1);

    // Factor x^m + 1
    Poly factors[128];
    int nfactors = irreducible_factors_xm_plus_1(m, factors);

    // DP over odd parts: map odd_lcm -> max two-exponent
    DPEntry dp[8192];
    int dp_count = 1;
    dp[0].odd_lcm = 1;
    dp[0].smax = 0;

    for (int fi = 0; fi < nfactors; fi++) {
        const Poly *p = &factors[fi];

        // g mod p
        Poly g_mod_p;
        poly_mod(&g_mod_p, &g_poly, p);

        // Check if g is invertible mod p
        Poly g_gcd;
        poly_gcd(&g_gcd, &g_mod_p, p);
        if (poly_deg(&g_gcd) > 0) continue;  // not invertible

        i64 odd_order = multiplicative_order(&g_mod_p, p);
        int smax = max_two_lift_exponent(&g_poly, odd_order, p, max_exp);

        // Update DP
        DPEntry new_dp[8192];
        int new_dp_count = 0;

        // Copy existing entries
        for (int i = 0; i < dp_count; i++) {
            new_dp[new_dp_count++] = dp[i];
        }

        // Add new entries by combining with this factor
        for (int i = 0; i < dp_count; i++) {
            i64 nl = ilcm(dp[i].odd_lcm, odd_order);
            int ns = dp[i].smax;
            if (smax > ns) ns = smax;

            // Check if (nl, ns) already in new_dp
            int found = 0;
            for (int j = 0; j < new_dp_count; j++) {
                if (new_dp[j].odd_lcm == nl) {
                    if (ns > new_dp[j].smax) new_dp[j].smax = ns;
                    found = 1;
                    break;
                }
            }
            if (!found) {
                new_dp[new_dp_count].odd_lcm = nl;
                new_dp[new_dp_count].smax = ns;
                new_dp_count++;
            }
        }

        memcpy(dp, new_dp, (size_t)new_dp_count * sizeof(DPEntry));
        dp_count = new_dp_count;
    }

    // Generate periods
    *out_count = 0;
    for (int i = 0; i < dp_count; i++) {
        for (int e = 0; e <= dp[i].smax; e++) {
            out_periods[(*out_count)++] = dp[i].odd_lcm << e;
        }
    }
}

// ==================== Main ====================

long long p843_native(void) {
    init_sieve();

    i64 all_periods[200000];
    int all_periods_count = 0;

    // Use a simple set for all periods (sorted insertion, no duplicates)
    periods_set_count = 0;

    for (int n = 3; n <= 100; n++) {
        i64 periods[10000];
        int np;
        periods_for_n(n, periods, &np);
        for (int i = 0; i < np; i++)
            periods_insert(periods[i]);
    }

    i64 total = 0;
    for (int i = 0; i < periods_set_count; i++)
        total += periods_set[i];

    return total;
}
