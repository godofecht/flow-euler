#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Project Euler 925: Next permutation of squares
   Sum of B(n^2) for n=1..10^16-1, mod 1e9+7.
   B(n) = next lexicographic permutation of n's digits (0 if none).

   Algorithm: square_sum_below_power10(k) + correction_sum(k).
   Port of the Python reference solver.

   Key challenge: numbers like 'visible' and 'representative' can be up to ~10^32,
   which fits in __int128 (max ~1.7*10^38). We compute delta_mod using i128. */

typedef int64_t i64;
typedef __int128 i128;

#define MOD 1000000007LL

/* ---------- delta_mod for i128 ---------- */
/* Compute (next_permutation(n) - n) % MOD for n given as i128 */
static i64 delta_mod_i128(i128 n) {
    if (n == 0) return 0;
    /* Extract digits */
    int digits[40];
    int len = 0;
    i128 tmp = n;
    while (tmp > 0) {
        digits[len++] = (int)(tmp % 10);
        tmp /= 10;
    }
    /* Reverse to get most-significant first */
    for (int i = 0, j = len - 1; i < j; i++, j--) {
        int t = digits[i]; digits[i] = digits[j]; digits[j] = t;
    }

    /* Compute n mod MOD from digits */
    i64 n_mod = 0;
    for (int d = 0; d < len; d++)
        n_mod = (i64)((i128)n_mod * 10 % MOD + digits[d]) % MOD;

    /* Find largest i with digits[i] < digits[i+1] */
    int i = len - 2;
    while (i >= 0 && digits[i] >= digits[i + 1]) i--;
    if (i < 0) {
        /* No next permutation, B(n) = 0, delta = -n mod MOD */
        i64 r = (MOD - n_mod) % MOD;
        return r;
    }

    /* Find largest j > i with digits[j] > digits[i] */
    int j = len - 1;
    while (digits[j] <= digits[i]) j--;

    /* Swap */
    int t = digits[i]; digits[i] = digits[j]; digits[j] = t;

    /* Reverse suffix starting at i+1 */
    int lo = i + 1, hi = len - 1;
    while (lo < hi) {
        t = digits[lo]; digits[lo] = digits[hi]; digits[hi] = t;
        lo++; hi--;
    }

    /* Compute next_permutation value mod MOD */
    i64 np_mod = 0;
    for (int d = 0; d < len; d++)
        np_mod = (i64)((i128)np_mod * 10 % MOD + digits[d]) % MOD;

    i64 dm = (np_mod - n_mod) % MOD;
    if (dm < 0) dm += MOD;
    return dm;
}

/* ---------- modular exponentiation ---------- */
static i64 mod_pow(i64 a, i64 e) {
    i64 r = 1; a %= MOD;
    if (a < 0) a += MOD;
    while (e > 0) {
        if (e & 1) r = (i64)((i128)r * a % MOD);
        a = (i64)((i128)a * a % MOD);
        e >>= 1;
    }
    return r;
}

/* ---------- square_sum_below_power10 ---------- */
static i64 square_sum_below_power10(int k) {
    i64 n = mod_pow(10, k);
    i64 t1 = (n - 1) % MOD; if (t1 < 0) t1 += MOD;
    i64 t2 = (2 * n - 1) % MOD; if (t2 < 0) t2 += MOD;
    i64 inv6 = mod_pow(6, MOD - 2);
    i64 result = (i64)((i128)n * t1 % MOD);
    result = (i64)((i128)result * t2 % MOD);
    result = (i64)((i128)result * inv6 % MOD);
    return result;
}

/* ---------- correction_sum ---------- */
/* Precomputed powers of 10 (exact, up to 2*k+2) and mod powers */
static i128 pow10[40];     /* exact powers, up to 2*k+2 */
static i64 pow10_mod[20];  /* mod powers, up to k */

static int K_global;

static i64 recurse(i64 suffix, int width, int trailing_zeros) {
    int k = K_global;

    if (width + trailing_zeros >= k) {
        /* delta_mod((suffix * pow10[trailing_zeros])^2) */
        i128 val = (i128)suffix * pow10[trailing_zeros];
        i128 sq = val * val;
        return delta_mod_i128(sq);
    }

    /* square_suffix = (suffix * suffix) % pow10[width] */
    i64 square_suffix = (i64)((i128)suffix * suffix % pow10[width]);
    i64 left_digit = square_suffix / (i64)pow10[width - 1];
    i128 next_place = pow10[width];
    i128 next_modulus = pow10[width + 1];
    int free_digits = k - width - trailing_zeros - 1;
    i64 completion_count = pow10_mod[free_digits];

    i64 total = 0;

    for (int digit = 0; digit < 10; digit++) {
        i128 next_suffix = (i128)digit * next_place + suffix;
        i128 next_sq = next_suffix * next_suffix;
        i64 new_digit = (i64)((next_sq % next_modulus) / next_place);

        if (new_digit >= left_digit) {
            total = (total + recurse((i64)next_suffix, width + 1, trailing_zeros)) % MOD;
        } else {
            /* visible = (next_square % next_modulus) * pow10[2*trailing_zeros] */
            i128 next_sq_mod = next_sq % next_modulus;
            i128 visible = next_sq_mod * pow10[2 * trailing_zeros];
            /* representative = pow10[width+1+2*trailing_zeros] + visible */
            i128 representative = pow10[width + 1 + 2 * trailing_zeros] + visible;
            i64 representative_delta = delta_mod_i128(representative);

            if (next_sq < next_modulus) {
                total = (total + delta_mod_i128(visible)) % MOD;
                i64 term = (completion_count - 1) % MOD; if (term < 0) term += MOD;
                term = (i64)((i128)term * representative_delta % MOD);
                total = (total + term) % MOD;
            } else {
                i64 term = (i64)((i128)completion_count * representative_delta % MOD);
                total = (total + term) % MOD;
            }
        }
    }
    return total;
}

static i64 correction_sum(int k) {
    K_global = k;

    /* Build pow10 array (exact, up to 2*k+2) using i128 */
    pow10[0] = 1;
    for (int i = 1; i <= 2 * k + 2; i++) pow10[i] = pow10[i - 1] * 10;

    /* Build pow10_mod (up to k) */
    pow10_mod[0] = 1;
    for (int i = 1; i <= k; i++) pow10_mod[i] = pow10_mod[i - 1] * 10 % MOD;

    i64 total = 0;
    for (int trailing_zeros = 0; trailing_zeros < k; trailing_zeros++) {
        for (int last_digit = 1; last_digit < 10; last_digit++) {
            total = (total + recurse(last_digit, 1, trailing_zeros)) % MOD;
        }
    }
    return total;
}

long long p925_native(void) {
    int k = 16;
    i64 result = (square_sum_below_power10(k) + correction_sum(k)) % MOD;
    return result;
}
