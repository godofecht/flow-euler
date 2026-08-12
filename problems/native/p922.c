#include <stdint.h>
#include <string.h>

typedef long long i64;

/* Project Euler 922 - Staircase Game
   R(m, w) mod 10^9+7, compute R(8, 64).
   Uses FWT for XOR convolution and polynomial exponentiation.
*/

#define MOD 1000000007LL
#define XOR_SIZE 64
#define W 64
#define M 8
#define DMAX (W - 2)         /* 62 */
#define DIFF_COUNT (2*DMAX+1) /* 125 */
#define MAX_POLY 1024

static i64 powmod(i64 a, i64 b, i64 m) {
    i64 r = 1; a %= m; if (a < 0) a += m;
    while (b > 0) {
        if (b & 1) r = (__int128)r * a % m;
        a = (__int128)a * a % m;
        b >>= 1;
    }
    return r;
}

static void fwt_xor(i64 *arr, int n, int inverse) {
    int step = 1;
    while (step < n) {
        int jump = step * 2;
        for (int i = 0; i < n; i += jump) {
            for (int j = i; j < i + step; j++) {
                i64 x = arr[j], y = arr[j + step];
                i64 s = x + y; if (s >= MOD) s -= MOD;
                i64 d = x - y; if (d < 0) d += MOD;
                arr[j] = s;
                arr[j + step] = d;
            }
        }
        step = jump;
    }
    if (inverse) {
        i64 inv_n = powmod(n, MOD - 2, MOD);
        for (int i = 0; i < n; i++)
            arr[i] = (__int128)arr[i] * inv_n % MOD;
    }
}

typedef struct { i64 coeffs[MAX_POLY]; int len; i64 offset; } poly_t;

static void poly_mul(poly_t *res, const poly_t *a, const poly_t *b) {
    const poly_t *big = a, *small = b;
    if (a->len < b->len) { big = b; small = a; }
    int new_len = big->len + small->len - 1;
    i64 tmp[MAX_POLY];
    memset(tmp, 0, new_len * sizeof(i64));
    for (int i = 0; i < big->len; i++) {
        if (big->coeffs[i] == 0) continue;
        for (int j = 0; j < small->len; j++) {
            if (small->coeffs[j] == 0) continue;
            tmp[i+j] = (tmp[i+j] + (__int128)big->coeffs[i] * small->coeffs[j]) % MOD;
        }
    }
    res->len = new_len;
    res->offset = a->offset + b->offset;
    memcpy(res->coeffs, tmp, new_len * sizeof(i64));
}

static void poly_pow(poly_t *res, const poly_t *base_in, int exp) {
    poly_t result;
    result.len = 1; result.offset = 0; result.coeffs[0] = 1;
    poly_t base = *base_in;
    poly_t tmp;
    while (exp > 0) {
        if (exp & 1) { poly_mul(&tmp, &result, &base); result = tmp; }
        exp >>= 1;
        if (exp) { poly_mul(&tmp, &base, &base); base = tmp; }
    }
    *res = result;
}

long long p922_native(void) {
    int w = W, m = M;
    int dmax = w - 2;
    int diff_count = 2 * dmax + 1;

    static i64 counts[DIFF_COUNT][XOR_SIZE];
    memset(counts, 0, sizeof(counts));

    for (int k = 1; k < w - 1; k++) {
        int limit = w - k;
        if (limit < 2) continue;
        int g = k - 1;
        int tmax = limit - 2;
        for (int t = 0; t <= tmax; t++) {
            int c = (limit - t) / 2;
            if (c <= 0) continue;
            counts[dmax + t][g] = (counts[dmax + t][g] + c) % MOD;
            if (t != 0)
                counts[dmax - t][g] = (counts[dmax - t][g] + c) % MOD;
        }
    }

    for (int d = 0; d < diff_count; d++)
        fwt_xor(counts[d], XOR_SIZE, 0);

    int final_offset = m * dmax;
    int final_len = 2 * final_offset + 1;

    static i64 Qhat[XOR_SIZE][MAX_POLY];
    memset(Qhat, 0, sizeof(Qhat));

    for (int t = 0; t < XOR_SIZE; t++) {
        poly_t poly;
        poly.len = diff_count;
        poly.offset = dmax;
        for (int d = 0; d < diff_count; d++)
            poly.coeffs[d] = counts[d][t];

        poly_t powered;
        poly_pow(&powered, &poly, m);

        for (int i = 0; i < powered.len && i < final_len; i++)
            Qhat[t][i] = powered.coeffs[i];
    }

    i64 ans = 0;
    i64 vec[XOR_SIZE];
    for (int idx = 0; idx < final_len; idx++) {
        for (int t = 0; t < XOR_SIZE; t++)
            vec[t] = Qhat[t][idx];
        fwt_xor(vec, XOR_SIZE, 1);

        int total_diff = idx - final_offset;
        if (total_diff > 0) {
            i64 s = 0;
            for (int t = 0; t < XOR_SIZE; t++) {
                s += vec[t]; if (s >= MOD) s -= MOD;
            }
            ans += s; if (ans >= MOD) ans -= MOD;
        } else if (total_diff == 0) {
            i64 s = 0;
            for (int t = 1; t < XOR_SIZE; t++) {
                s += vec[t]; if (s >= MOD) s -= MOD;
            }
            ans += s; if (ans >= MOD) ans -= MOD;
        }
    }

    return ans;
}
