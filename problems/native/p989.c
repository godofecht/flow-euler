// Project Euler 989
// Fibonacci sum with nonprimitive pairs, mod 1e9+9.
// Port of the Python reference solver to C.
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef int64_t i64;
typedef __int128 i128;

#define MOD 1000000009LL
#define TARGET_LIMIT 100000000000000LL  /* 10^14 */
#define SMALL_LIMIT 8

static i64 mulmod(i64 a, i64 b) {
    return (i64)((i128)a * b % MOD);
}

static i64 mod_pow(i64 a, i64 e) {
    i64 r = 1;
    a %= MOD;
    if (a < 0) a += MOD;
    while (e > 0) {
        if (e & 1) r = mulmod(r, a);
        a = mulmod(a, a);
        e >>= 1;
    }
    return r;
}

static i64 tonelli_shanks(i64 n, i64 p) {
    if (n == 0) return 0;
    if (mod_pow(n, (p - 1) / 2) != 1) return -1;
    if (p % 4 == 3) return mod_pow(n, (p + 1) / 4);

    i64 q = p - 1;
    i64 s = 0;
    while (q % 2 == 0) { q /= 2; s++; }

    i64 z = 2;
    while (mod_pow(z, (p - 1) / 2) != p - 1) z++;

    i64 m = s;
    i64 c = mod_pow(z, q);
    i64 t = mod_pow(n, q);
    i64 r = mod_pow(n, (q + 1) / 2);

    while (t != 1) {
        i64 i = 1;
        i64 t2i = mulmod(t, t);
        while (t2i != 1) { t2i = mulmod(t2i, t2i); i++; }
        i64 b = mod_pow(c, (i64)1 << (m - i - 1));
        r = mulmod(r, b);
        c = mulmod(b, b);
        t = mulmod(t, c);
        m = i;
    }
    return r;
}

static i64 isqrt_i64(i64 n) {
    if (n <= 0) return 0;
    i64 r = (i64)sqrt((double)n);
    while (r > 0 && r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r;
}

static int8_t *mobius_sieve(i64 limit) {
    int8_t *mu = (int8_t *)malloc((size_t)(limit + 1));
    memset(mu, 1, (size_t)(limit + 1));
    uint8_t *is_prime = (uint8_t *)malloc((size_t)(limit + 1));
    memset(is_prime, 1, (size_t)(limit + 1));
    if (limit >= 0) is_prime[0] = 0;
    if (limit >= 1) is_prime[1] = 0;

    for (i64 p = 2; p <= limit; p++) {
        if (!is_prime[p]) continue;
        for (i64 mult = p; mult <= limit; mult += p)
            mu[mult] = -mu[mult];
        i64 square = p * p;
        if (square <= limit) {
            for (i64 mult = square; mult <= limit; mult += square)
                mu[mult] = 0;
            for (i64 mult = square; mult <= limit; mult += p)
                is_prime[mult] = 0;
        }
        for (i64 mult = p + p; mult <= limit; mult += p)
            is_prime[mult] = 0;
    }
    free(is_prime);
    return mu;
}

/* Small nonprimitive terms: q = a^2 - ab - b^2, 0 < q <= 8 */
static int small_values[16];
static int small_value_count = 0;

static void build_small_values(void) {
    int max_a = 2 * (int)isqrt_i64(SMALL_LIMIT) + 2;
    for (int a = 2; a <= max_a; a++) {
        for (int b = 1; b <= a / 2; b++) {
            int q = a * a - a * b - b * b;
            if (q > 0 && q <= SMALL_LIMIT)
                small_values[small_value_count++] = q;
        }
    }
    for (int i = 0; i < small_value_count; i++)
        for (int j = i + 1; j < small_value_count; j++)
            if (small_values[j] < small_values[i]) {
                int t = small_values[i];
                small_values[i] = small_values[j];
                small_values[j] = t;
            }
}

static void eval_small_nonprimitive_pair(i64 limit, i64 z1, i64 z2,
                                          i64 *out1, i64 *out2) {
    i64 total1 = 0, total2 = 0;
    i64 power1 = 1, power2 = 1;
    i64 exponent = 0;

    for (int idx = 0; idx < small_value_count; idx++) {
        i64 target = small_values[idx];
        if (target > limit) break;
        while (exponent < target) {
            power1 = mulmod(power1, z1);
            power2 = mulmod(power2, z2);
            exponent++;
        }
        total1 += power1;
        if (total1 >= MOD) total1 -= MOD;
        total2 += power2;
        if (total2 >= MOD) total2 -= MOD;
    }
    *out1 = total1;
    *out2 = total2;
}

static void nonprimitive_pair(i64 limit, i64 z1, i64 z1_inv,
                               i64 z2, i64 z2_inv,
                               i64 *out1, i64 *out2) {
    if (limit <= SMALL_LIMIT) {
        eval_small_nonprimitive_pair(limit, z1, z2, out1, out2);
        return;
    }

    i64 total1 = 0, total2 = 0;

    i64 z1_sq = mulmod(z1, z1);
    i64 z2_sq = mulmod(z2, z2);

    i64 z1_inv_sq = mulmod(z1_inv, z1_inv);
    i64 z2_inv_sq = mulmod(z2_inv, z2_inv);
    i64 z1_inv_4 = mulmod(z1_inv_sq, z1_inv_sq);
    i64 z2_inv_4 = mulmod(z2_inv_sq, z2_inv_sq);
    i64 z1_inv_5 = mulmod(z1_inv_4, z1_inv);
    i64 z2_inv_5 = mulmod(z2_inv_4, z2_inv);
    i64 z1_inv_10 = mulmod(z1_inv_5, z1_inv_5);
    i64 z2_inv_10 = mulmod(z2_inv_5, z2_inv_5);
    i64 z1_inv_15 = mulmod(z1_inv_10, z1_inv_5);
    i64 z2_inv_15 = mulmod(z2_inv_10, z2_inv_5);

    /* Even part: sliding window of z^(m^2) */
    i64 even_weight1 = z1_inv_5, even_weight2 = z2_inv_5;
    i64 even_delta1 = z1_inv_15, even_delta2 = z2_inv_15;

    i64 add_index = 0;
    i64 add_term1 = 1, add_term2 = 1;
    i64 add_step1 = z1, add_step2 = z2;

    i64 drop_index = 0;
    i64 drop_term1 = 1, drop_term2 = 1;
    i64 drop_step1 = z1, drop_step2 = z2;

    i64 window1 = 0, window2 = 0;
    i64 t = 1;
    i64 lower = 3, upper = 0;
    i64 rhs = limit + 5;

    while ((upper + 1) * (upper + 1) <= rhs) upper++;

    while (lower <= upper) {
        while (add_index <= upper) {
            window1 += add_term1; if (window1 >= MOD) window1 -= MOD;
            window2 += add_term2; if (window2 >= MOD) window2 -= MOD;
            add_term1 = mulmod(add_term1, add_step1);
            add_step1 = mulmod(add_step1, z1_sq);
            add_term2 = mulmod(add_term2, add_step2);
            add_step2 = mulmod(add_step2, z2_sq);
            add_index++;
        }
        while (drop_index < lower) {
            window1 -= drop_term1; if (window1 < 0) window1 += MOD;
            window2 -= drop_term2; if (window2 < 0) window2 += MOD;
            drop_term1 = mulmod(drop_term1, drop_step1);
            drop_step1 = mulmod(drop_step1, z1_sq);
            drop_term2 = mulmod(drop_term2, drop_step2);
            drop_step2 = mulmod(drop_step2, z2_sq);
            drop_index++;
        }
        total1 = (total1 + mulmod(window1, even_weight1)) % MOD;
        total2 = (total2 + mulmod(window2, even_weight2)) % MOD;
        even_weight1 = mulmod(even_weight1, even_delta1);
        even_delta1 = mulmod(even_delta1, z1_inv_10);
        even_weight2 = mulmod(even_weight2, even_delta2);
        even_delta2 = mulmod(even_delta2, z2_inv_10);
        rhs += 10 * t + 5;
        t++;
        lower += 3;
        while ((upper + 1) * (upper + 1) <= rhs) upper++;
    }

    /* Odd part: sliding window of z^(m(m+1)) */
    i64 odd_weight1 = z1_inv, odd_weight2 = z2_inv;
    i64 odd_delta1 = z1_inv_10, odd_delta2 = z2_inv_10;

    add_index = 0;
    add_term1 = 1; add_term2 = 1;
    add_step1 = z1_sq; add_step2 = z2_sq;
    drop_index = 0;
    drop_term1 = 1; drop_term2 = 1;
    drop_step1 = z1_sq; drop_step2 = z2_sq;
    window1 = 0; window2 = 0;
    t = 0;
    lower = 1; upper = 0;
    rhs = limit + 1;

    while ((upper + 1) * (upper + 2) <= rhs) upper++;

    while (lower <= upper) {
        while (add_index <= upper) {
            window1 += add_term1; if (window1 >= MOD) window1 -= MOD;
            window2 += add_term2; if (window2 >= MOD) window2 -= MOD;
            add_term1 = mulmod(add_term1, add_step1);
            add_step1 = mulmod(add_step1, z1_sq);
            add_term2 = mulmod(add_term2, add_step2);
            add_step2 = mulmod(add_step2, z2_sq);
            add_index++;
        }
        while (drop_index < lower) {
            window1 -= drop_term1; if (window1 < 0) window1 += MOD;
            window2 -= drop_term2; if (window2 < 0) window2 += MOD;
            drop_term1 = mulmod(drop_term1, drop_step1);
            drop_step1 = mulmod(drop_step1, z1_sq);
            drop_term2 = mulmod(drop_term2, drop_step2);
            drop_step2 = mulmod(drop_step2, z2_sq);
            drop_index++;
        }
        total1 = (total1 + mulmod(window1, odd_weight1)) % MOD;
        total2 = (total2 + mulmod(window2, odd_weight2)) % MOD;
        odd_weight1 = mulmod(odd_weight1, odd_delta1);
        odd_delta1 = mulmod(odd_delta1, z1_inv_10);
        odd_weight2 = mulmod(odd_weight2, odd_delta2);
        odd_delta2 = mulmod(odd_delta2, z2_inv_10);
        rhs += 10 * t + 10;
        t++;
        lower += 3;
        while ((upper + 1) * (upper + 2) <= rhs) upper++;
    }

    *out1 = total1;
    *out2 = total2;
}

long long p989_native(void) {
    build_small_values();

    i64 limit = TARGET_LIMIT;
    i64 root = isqrt_i64(limit);
    int8_t *mu = mobius_sieve(root);

    i64 sqrt5_mod = tonelli_shanks(5, MOD);
    i64 inv_sqrt5_mod = mod_pow(sqrt5_mod, MOD - 2);
    i64 inv2_mod = (MOD + 1) / 2;
    i64 phi_mod = mulmod((1 + sqrt5_mod) % MOD, inv2_mod);
    i64 phi_inv_mod = mod_pow(phi_mod, MOD - 2);
    i64 phi_sq_mod = mulmod(phi_mod, phi_mod);
    i64 phi_inv_sq_mod = mulmod(phi_inv_mod, phi_inv_mod);

    i64 p_phi = 0, p_psi = 0;
    i64 phi_pow_g2 = 1, phi_inv_pow_g2 = 1;
    i64 forward_step = phi_mod, backward_step = phi_inv_mod;
    i64 g_square = 1;

    for (i64 g = 1; g <= root; g++) {
        phi_pow_g2 = mulmod(phi_pow_g2, forward_step);
        forward_step = mulmod(forward_step, phi_sq_mod);
        phi_inv_pow_g2 = mulmod(phi_inv_pow_g2, backward_step);
        backward_step = mulmod(backward_step, phi_inv_sq_mod);

        int8_t mu_g = mu[g];
        if (mu_g) {
            i64 scaled_limit = limit / g_square;
            i64 psi_pow_g2, psi_inv_pow_g2;
            if (g & 1) {
                psi_pow_g2 = MOD - phi_inv_pow_g2;
                psi_inv_pow_g2 = MOD - phi_pow_g2;
            } else {
                psi_pow_g2 = phi_inv_pow_g2;
                psi_inv_pow_g2 = phi_pow_g2;
            }

            i64 np_phi, np_psi;
            nonprimitive_pair(scaled_limit,
                              phi_pow_g2, phi_inv_pow_g2,
                              psi_pow_g2, psi_inv_pow_g2,
                              &np_phi, &np_psi);

            if (mu_g == 1) {
                p_phi += np_phi; if (p_phi >= MOD) p_phi -= MOD;
                p_psi += np_psi; if (p_psi >= MOD) p_psi -= MOD;
            } else {
                p_phi -= np_phi; if (p_phi < 0) p_phi += MOD;
                p_psi -= np_psi; if (p_psi < 0) p_psi += MOD;
            }
        }
        g_square += 2 * g + 1;
    }

    free(mu);
    return mulmod((p_phi - p_psi + MOD) % MOD, inv_sqrt5_mod);
}
