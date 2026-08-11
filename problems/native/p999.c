/* Project Euler 999
 * Rescaled elliptic divisibility sequence: a_n mod 1234567891.
 * Answer = a_{10^18 + 3} = 801096743.
 *
 * Uses the doubling recurrence for the EDS block (W_{n-3}..W_{n+4})
 * in O(log n), then rescales by INV_TWO^(n^2/4).
 */
#include <stdint.h>
#include <stdio.h>

typedef int64_t i64;
typedef __int128 i128;

#define MOD 1234567891LL
#define INV_TWO 617283946LL  /* (MOD + 1) / 2, inverse of 2 mod MOD */

static const i64 SMALL_W[9] = {0, 1, 2, -4, -32, -192, 3584, 77824, 262144};

static i64 norm(i64 x) {
    x %= MOD;
    if (x < 0) x += MOD;
    return x;
}

static i64 mod_mul(i64 a, i64 b) {
    return (norm(a) * norm(b)) % MOD;  /* both < 1.23e9, product < 1.5e18 */
}

static i64 mod_sub(i64 a, i64 b) {
    i64 r = norm(a) - norm(b);
    if (r < 0) r += MOD;
    return r;
}

static i64 mod_pow_i64(i64 base, i64 exp) {
    i64 b = norm(base), result = 1;
    while (exp > 0) {
        if (exp & 1) result = mod_mul(result, b);
        b = mod_mul(b, b);
        exp >>= 1;
    }
    return result;
}

static i64 mod_pow_i128(i64 base, i128 exp) {
    i64 b = norm(base), result = 1;
    while (exp > 0) {
        if (exp & 1) result = mod_mul(result, b);
        b = mod_mul(b, b);
        exp >>= 1;
    }
    return result;
}

static i64 small_w(i64 index) {
    if (index < 0) return norm(-small_w(-index));
    return norm(SMALL_W[index]);
}

typedef struct { i64 v[8]; } Block;

/* get(index) into the source block whose first element is W_{source_start}. */
static i64 blk_get(Block s, i64 start, i64 index) {
    return s.v[index - start];
}

/* W_{2*index - 1} from terms around W_index. */
static i64 odd_val(Block s, i64 start, i64 index) {
    i64 a = blk_get(s, start, index + 1);
    i64 b = mod_pow_i64(blk_get(s, start, index - 1), 3);
    i64 c = blk_get(s, start, index - 2);
    i64 d = mod_pow_i64(blk_get(s, start, index), 3);
    return mod_sub(mod_mul(a, b), mod_mul(c, d));
}

/* W_{2*index}; division is only by the fixed W_2 = 2. */
static i64 even_val(Block s, i64 start, i64 index) {
    i64 a = blk_get(s, start, index);
    i64 p1 = mod_pow_i64(blk_get(s, start, index - 1), 2);
    i64 c2 = blk_get(s, start, index + 2);
    i64 p2 = mod_pow_i64(blk_get(s, start, index + 1), 2);
    i64 e2 = blk_get(s, start, index - 2);
    i64 inner = mod_sub(mod_mul(c2, p1), mod_mul(e2, p2));
    return mod_mul(mod_mul(a, INV_TWO), inner);
}

/* Return (W_{n-3}, ..., W_{n+4}) mod MOD. */
static Block eds_block(i64 n) {
    Block r;
    if (n <= 4) {
        for (int i = 0; i < 8; i++) r.v[i] = small_w(n - 3 + i);
        return r;
    }
    i64 middle = n / 2;
    Block s = eds_block(middle);
    i64 start = middle - 3;
    if (n % 2 == 0) {
        r.v[0] = odd_val(s, start, middle - 1);
        r.v[1] = even_val(s, start, middle - 1);
        r.v[2] = odd_val(s, start, middle);
        r.v[3] = even_val(s, start, middle);
        r.v[4] = odd_val(s, start, middle + 1);
        r.v[5] = even_val(s, start, middle + 1);
        r.v[6] = odd_val(s, start, middle + 2);
        r.v[7] = even_val(s, start, middle + 2);
    } else {
        r.v[0] = even_val(s, start, middle - 1);
        r.v[1] = odd_val(s, start, middle);
        r.v[2] = even_val(s, start, middle);
        r.v[3] = odd_val(s, start, middle + 1);
        r.v[4] = even_val(s, start, middle + 1);
        r.v[5] = odd_val(s, start, middle + 2);
        r.v[6] = even_val(s, start, middle + 2);
        r.v[7] = odd_val(s, start, middle + 3);
    }
    return r;
}

static i64 w_mod(i64 n) {
    if (n < 0) return norm(-w_mod(-n));
    return eds_block(n).v[3];
}

static i64 a_mod(i64 n) {
    i64 rem = n % 4;
    i64 sign = (rem == 1 || rem == 2) ? 1 : -1;
    i128 nn = (i128)n * (i128)n;
    i128 exp = nn / 4;
    i64 inverse_scale = mod_pow_i128(INV_TWO, exp);
    i64 w = w_mod(n);
    i64 sw = norm(sign * w);
    return mod_mul(sw, inverse_scale);
}

long long p999_native(void) {
    return a_mod((i64)1000000000000000003LL);  /* 10^18 + 3 */
}
