/*
 * Project Euler 792 — "Too Many Twos"
 *
 * S(n) = sum_{k=1..n} (-2)^k * C(2k,k)
 * u(n) = v2(3*S(n) + 4)          (2-adic valuation)
 * U(N) = sum_{n=1..N} u(n^3),    N = 10000.
 *
 * In the 2-adics 3*sum_{k>=1} (-2)^k*C(2k,k) + 4 = 0, so
 *   3*S(n) + 4 = -3 * sum_{k>n} R(k),  R(k) = (-2)^k * C(2k,k).
 * -3 is odd, so u(n) = v2( sum_{k>n} R(k) ).
 *
 * R(k+1)/R(k) = -4 * (2k+1)/(k+1),  v2(R(k)) = k + popcount(k).
 *
 * We sum the tail k = n+1 .. n+m (m up to 220) keeping only odd parts
 * modulo 2^P (powers of two tracked separately).  The remainder for
 * k >= n+m+1 is divisible by 2^(n+m+1), so once the partial sum's
 * valuation drops below n+m+1 it equals u(n).  Precision P starts at
 * 256 bits and is doubled when the reduced sum vanishes mod 2^P.
 */

#include <stdint.h>
#include <stdio.h>
#include <string.h>

typedef uint64_t u64;
typedef int64_t  i64;
typedef unsigned __int128 u128;

#define MAXLIMBS 32   /* 2048 bits */

/* ---- bignum mod 2^(64*nl), little-endian limbs ---- */

static int bn_is_zero(const u64 *a, int nl) {
    for (int i = 0; i < nl; i++) if (a[i]) return 0;
    return 1;
}

/* trailing zero bits; returns nl*64 if all zero */
static int bn_v2(const u64 *a, int nl) {
    for (int i = 0; i < nl; i++)
        if (a[i]) return i * 64 + __builtin_ctzll(a[i]);
    return nl * 64;
}

static void bn_zero(u64 *a, int nl) { memset(a, 0, (size_t)nl * sizeof(u64)); }

static void bn_copy(u64 *dst, const u64 *src, int nl) {
    memmove(dst, src, (size_t)nl * sizeof(u64));
}

static void bn_from_u64(u64 *a, u64 v, int nl) {
    bn_zero(a, nl);
    a[0] = v;
}

/* a = (i64)v sign-extended mod 2^(64*nl) */
static void bn_from_i64(u64 *a, i64 v, int nl) {
    u64 fill = (v < 0) ? ~(u64)0 : (u64)0;
    for (int i = 1; i < nl; i++) a[i] = fill;
    a[0] = (u64)v;
}

/* dst = (a + b) mod 2^(64*nl) */
static void bn_add(u64 *dst, const u64 *a, const u64 *b, int nl) {
    u64 carry = 0;
    for (int i = 0; i < nl; i++) {
        u128 s = (u128)a[i] + b[i] + carry;
        dst[i] = (u64)s;
        carry = (u64)(s >> 64);
    }
}

/* dst = (a << shift) mod 2^(64*nl) */
static void bn_shl(u64 *dst, const u64 *a, int shift, int nl) {
    if (shift <= 0) { bn_copy(dst, a, nl); return; }
    if (shift >= nl * 64) { bn_zero(dst, nl); return; }
    int ws = shift / 64, bs = shift % 64;
    bn_zero(dst, nl);
    if (bs == 0) {
        for (int i = nl - 1; i >= ws; i--) dst[i] = a[i - ws];
    } else {
        for (int i = nl - 1; i >= ws; i--) {
            u64 lo = a[i - ws] << bs;
            u64 hi = (i - ws - 1 >= 0) ? (a[i - ws - 1] >> (64 - bs)) : 0;
            dst[i] = lo | hi;
        }
    }
}

/* dst = (a * b) mod 2^(64*nl) */
static void bn_mul(u64 *dst, const u64 *a, const u64 *b, int nl) {
    u64 tmp[MAXLIMBS];
    bn_zero(tmp, nl);
    for (int i = 0; i < nl; i++) {
        u64 carry = 0;
        for (int j = 0; j < nl - i; j++) {
            u128 p = (u128)a[i] * b[j] + tmp[i + j] + carry;
            tmp[i + j] = (u64)p;
            carry = (u64)(p >> 64);
        }
    }
    bn_copy(dst, tmp, nl);
}

/* modular inverse of odd a (low limb set, rest zero) mod 2^(64*nl), Newton */
static void bn_inv_odd(u64 *dst, const u64 *a, int nl) {
    u64 x[MAXLIMBS], t[MAXLIMBS], two[MAXLIMBS];
    bn_from_u64(x, 1, nl);          /* correct mod 2 for any odd a */
    bn_from_u64(two, 2, nl);
    int bits = nl * 64;
    for (int it = 0; (1 << it) < bits; it++) {
        /* x = x * (2 - a*x) */
        u64 ax[MAXLIMBS];
        bn_mul(ax, a, x, nl);
        u64 s[MAXLIMBS];
        /* s = 2 - ax  (mod 2^(64*nl)) = (~ax + 1) + 2 ... use: two - ax */
        bn_zero(s, nl);
        /* s = two - ax mod 2^(64*nl): compute ax_neg = ~ax + 1, then + two */
        u64 neg[MAXLIMBS];
        for (int i = 0; i < nl; i++) neg[i] = ~ax[i];
        /* add 1 */
        u64 c = 1;
        for (int i = 0; i < nl && c; i++) {
            u128 s2 = (u128)neg[i] + c;
            neg[i] = (u64)s2;
            c = (u64)(s2 >> 64);
        }
        bn_add(s, neg, two, nl);
        bn_mul(t, x, s, nl);
        bn_copy(x, t, nl);
    }
    bn_copy(dst, x, nl);
}

/* ---- per-call inverse cache (denom_odd fits in u64) ---- */

struct inv_entry { u64 key; u64 val[MAXLIMBS]; };

static struct inv_entry inv_cache[256];
static int inv_cache_n;

static const u64 *inv_odd(u64 denom_odd, int nl) {
    for (int i = 0; i < inv_cache_n; i++)
        if (inv_cache[i].key == denom_odd) return inv_cache[i].val;
    /* compute and store */
    u64 a[MAXLIMBS];
    bn_from_u64(a, denom_odd, nl);
    bn_inv_odd(inv_cache[inv_cache_n].val, a, nl);
    inv_cache[inv_cache_n].key = denom_odd;
    return inv_cache[inv_cache_n++].val;
}

/* ---- core: u(n) ---- */

static i64 u_of(i64 n) {
    static const int P_list[] = {256, 512, 1024, 2048};
    for (int pi = 0; pi < 4; pi++) {
        int P = P_list[pi];
        int nl = P / 64;
        inv_cache_n = 0;

        i64 k = n + 1;
        i64 exp = k + (i64)__builtin_popcountll((u64)k);

        u64 odd[MAXLIMBS];
        bn_from_u64(odd, 1, nl);

        int have_min = 0;
        i64 min_exp = 0;
        u64 scaled_sum[MAXLIMBS];
        bn_zero(scaled_sum, nl);

        for (int m = 1; m <= 220; m++) {
            if (!have_min) {
                min_exp = exp;
                bn_copy(scaled_sum, odd, nl);
                have_min = 1;
            } else if (exp < min_exp) {
                i64 shift = min_exp - exp;
                if (shift < P) bn_shl(scaled_sum, scaled_sum, (int)shift, nl);
                else           bn_zero(scaled_sum, nl);
                min_exp = exp;
                bn_add(scaled_sum, scaled_sum, odd, nl);
            } else {
                i64 shift = exp - min_exp;
                if (shift < P) {
                    u64 tmp[MAXLIMBS];
                    bn_shl(tmp, odd, (int)shift, nl);
                    bn_add(scaled_sum, scaled_sum, tmp, nl);
                }
            }

            if (bn_is_zero(scaled_sum, nl)) break;  /* need more precision */

            i64 v_partial = min_exp + (i64)bn_v2(scaled_sum, nl);
            if (v_partial < n + m + 1) return v_partial;

            /* advance R(k) -> R(k+1): ratio = -4*(2k+1)/(k+1) */
            i64 denom = k + 1;
            int t = __builtin_ctzll((u64)denom);
            u64 denom_odd = (u64)denom >> t;

            const u64 *inv = inv_odd(denom_odd, nl);
            u64 factor[MAXLIMBS];
            bn_from_i64(factor, -(2 * k + 1), nl);  /* -(2k+1) mod 2^P */
            u64 tmp2[MAXLIMBS];
            bn_mul(tmp2, factor, inv, nl);
            u64 odd_new[MAXLIMBS];
            bn_mul(odd_new, odd, tmp2, nl);
            bn_copy(odd, odd_new, nl);

            exp += 2 - t;
            k = denom;

            /* consistency: exp == k + popcount(k) */
            if (exp != k + (i64)__builtin_popcountll((u64)k)) {
                /* internal error; should never happen */
                return -1;
            }
        }
        /* try higher P */
    }
    return -1;  /* failed */
}

long long p792_native(void) {
    i64 total = 0;
    for (i64 n = 1; n <= 10000; n++) {
        i64 nc = n * n * n;
        total += u_of(nc);
    }
    return (long long)total;
}
