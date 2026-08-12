// Project Euler 996: Overtaking
// Polynomial generating function approach, mod 1234567891.
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef long long i64;
typedef __int128 i128;

#define MOD 1234567891LL

static i64 mod_add(i64 a, i64 b) {
    a += b;
    if (a >= MOD) a -= MOD;
    if (a < 0) a += MOD;
    return a;
}

static i64 mod_sub(i64 a, i64 b) {
    a -= b;
    if (a < 0) a += MOD;
    return a;
}

static i64 mod_mul(i64 a, i64 b) {
    return (i64)((i128)a * b % MOD);
}

static i64 ext_gcd(i64 a, i64 b, i64 *x, i64 *y) {
    if (b == 0) { *x = 1; *y = 0; return a; }
    i64 x1, y1;
    i64 g = ext_gcd(b, a % b, &x1, &y1);
    *x = y1;
    *y = x1 - (a / b) * y1;
    return g;
}

static i64 mod_inv(i64 a) {
    a %= MOD;
    if (a < 0) a += MOD;
    i64 x, y;
    ext_gcd(a, MOD, &x, &y);
    return ((x % MOD) + MOD) % MOD;
}

// Compute C(n, k) mod MOD for small k.
static i64 comb_mod(i64 n, i64 k) {
    if (k < 0 || k > n || n < 0) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;
    i64 num = 1, den = 1;
    for (i64 i = 0; i < k; i++) {
        num = mod_mul(num, (n - i) % MOD);
        den = mod_mul(den, (i + 1) % MOD);
    }
    return mod_mul(num, mod_inv(den));
}

// Polynomial operations. Polynomials are arrays of i64 coefficients.
// We use dynamic arrays with explicit lengths.

static void trim(i64 *poly, int *len) {
    while (*len > 1 && poly[*len - 1] == 0) (*len)--;
}

static void add_to(i64 *dst, int *dst_len, const i64 *src, int src_len) {
    if (*dst_len < src_len) {
        for (int i = *dst_len; i < src_len; i++) dst[i] = 0;
        *dst_len = src_len;
    }
    for (int i = 0; i < src_len; i++)
        dst[i] = mod_add(dst[i], src[i]);
}

// Multiply poly by (1 - q), return result in out.
static void mul_one_minus_q(const i64 *poly, int plen, i64 *out, int *out_len) {
    int n = plen + 1;
    for (int i = 0; i < n; i++) out[i] = 0;
    for (int i = 0; i < plen; i++) {
        out[i] = mod_add(out[i], poly[i]);
        out[i + 1] = mod_sub(out[i + 1], poly[i]);
    }
    *out_len = n;
    trim(out, out_len);
}

// Multiply a and b, discarding terms above max_degree.
static void mul_poly(const i64 *a, int alen, const i64 *b, int blen,
                     int max_degree, i64 *out, int *out_len) {
    if (alen == 0 || blen == 0) { out[0] = 0; *out_len = 1; return; }
    int n = (alen + blen - 2 < max_degree ? alen + blen - 2 : max_degree) + 1;
    for (int i = 0; i < n; i++) out[i] = 0;
    for (int i = 0; i < alen; i++) {
        if (a[i] == 0) continue;
        int last_j = (blen - 1 < max_degree - i ? blen - 1 : max_degree - i);
        for (int j = 0; j <= last_j; j++) {
            if (b[j])
                out[i + j] = mod_add(out[i + j], mod_mul(a[i], b[j]));
        }
    }
    *out_len = n;
    trim(out, out_len);
}

static i64 block_count(int length, int cost) {
    if (cost <= 0 || 2 * cost < length) return 0;
    i64 total = comb_mod(2 * cost - 1, length - 1);
    i64 too_large = (cost < length) ? 0 : comb_mod(cost - 1, length - 1);
    return mod_sub(total, mod_mul(length, too_large));
}

static void block_numerator(int length, i64 *out, int *out_len) {
    for (int j = 0; j <= length; j++) {
        i64 value = 0;
        for (int i = 0; i <= j; i++) {
            i64 sign = (i % 2) ? -1 : 1;
            i64 term = mod_mul(comb_mod(length, i), block_count(length, j - i));
            if (sign < 0) value = mod_sub(value, term);
            else value = mod_add(value, term);
        }
        out[j] = value;
    }
    *out_len = length + 1;
    trim(out, out_len);
}

// Main DP: compute numerator polynomial for all valid vectors of length n.
static void numerator_for_all_valid_vectors(int n, i64 *result, int *result_len) {
    // block_num[length] for length 2..n
    i64 *block_num[128];
    int block_num_len[128];
    for (int len = 2; len <= n; len++) {
        block_num[len] = calloc(len + 2, sizeof(i64));
        block_numerator(len, block_num[len], &block_num_len[len]);
    }

    // total[pos] and zero_end[pos]: arrays of polynomials
    i64 *total[128];
    int total_len[128];
    i64 *zero_end[128];
    int zero_end_len[128];

    for (int pos = 0; pos <= n; pos++) {
        total[pos] = calloc(n + 2, sizeof(i64));
        total_len[pos] = 1;
        total[pos][0] = 0;
        zero_end[pos] = calloc(n + 2, sizeof(i64));
        zero_end_len[pos] = 1;
        zero_end[pos][0] = 0;
    }
    total[0][0] = 1; total_len[0] = 1;
    zero_end[0][0] = 1; zero_end_len[0] = 1;

    i64 *tmp = calloc(n + 2, sizeof(i64));
    i64 *prod = calloc(n + 2, sizeof(i64));
    int tmp_len, prod_len;

    for (int pos = 0; pos <= n; pos++) {
        if (pos < n && total_len[pos] > 0 && !(total_len[pos] == 1 && total[pos][0] == 0)) {
            mul_one_minus_q(total[pos], total_len[pos], tmp, &tmp_len);
            add_to(total[pos + 1], &total_len[pos + 1], tmp, tmp_len);
            add_to(zero_end[pos + 1], &zero_end_len[pos + 1], tmp, tmp_len);
        }

        if (zero_end_len[pos] > 0 && !(zero_end_len[pos] == 1 && zero_end[pos][0] == 0)) {
            for (int length = 2; length <= n - pos; length++) {
                mul_poly(zero_end[pos], zero_end_len[pos],
                         block_num[length], block_num_len[length],
                         pos + length, prod, &prod_len);
                add_to(total[pos + length], &total_len[pos + length], prod, prod_len);
            }
        }
    }

    memcpy(result, total[n], total_len[n] * sizeof(i64));
    *result_len = total_len[n];

    for (int len = 2; len <= n; len++) free(block_num[len]);
    for (int pos = 0; pos <= n; pos++) { free(total[pos]); free(zero_end[pos]); }
    free(tmp); free(prod);
}

static i64 count_tuples(int n, i64 k) {
    i64 max_cost = k / 2;
    i64 *numer = calloc(n + 2, sizeof(i64));
    int numer_len;
    numerator_for_all_valid_vectors(n, numer, &numer_len);

    i64 answer = 0;
    for (int degree = 0; degree < numer_len; degree++) {
        if (numer[degree] == 0 || degree > max_cost) continue;
        i64 ways = comb_mod(max_cost - degree + n, n);
        answer = mod_add(answer, mod_mul(numer[degree], ways));
    }
    free(numer);
    return answer;
}

long long p996_native(void) {
    return count_tuples(123, 4567891);
}
