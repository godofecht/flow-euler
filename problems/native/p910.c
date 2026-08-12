// Project Euler 910: L-Expressions II
// Polynomial arithmetic in Z/(10^9)[x] mod prod_{k=1}^{40}(x+k).
// All polynomials are degree-39 remainders.
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define DEGREE 40
#define MOD 1000000000LL

static long long FOLD[DEGREE];
static long long ZERO[DEGREE];
static long long ONE[DEGREE];
static long long X[DEGREE];

static void init_globals(void) {
    // Build FOLD: coefficients to replace x^40
    // prod_{k=1}^{40} (x + k) = x^40 + c39*x^39 + ... + c0
    // x^40 = -(c39*x^39 + ... + c0) mod the product
    long long coeffs[DEGREE + 1];
    memset(coeffs, 0, sizeof(coeffs));
    coeffs[0] = 1;
    int len = 1;
    for (int k = 1; k <= DEGREE; k++) {
        long long nxt[DEGREE + 2];
        memset(nxt, 0, sizeof(nxt));
        for (int i = 0; i < len; i++) {
            nxt[i] = (nxt[i] + coeffs[i] * k) % MOD;
            nxt[i + 1] = (nxt[i + 1] + coeffs[i]) % MOD;
        }
        len++;
        memcpy(coeffs, nxt, len * sizeof(long long));
    }
    // FOLD[i] = -coeffs[i] for i=0..39
    for (int i = 0; i < DEGREE; i++) {
        FOLD[i] = ((MOD - coeffs[i]) % MOD);
    }

    memset(ZERO, 0, sizeof(ZERO));
    memset(ONE, 0, sizeof(ONE));
    ONE[0] = 1;
    memset(X, 0, sizeof(X));
    X[1] = 1;
}

static void poly_add(long long *a, long long *b, long long *out) {
    for (int i = 0; i < DEGREE; i++) {
        out[i] = (a[i] + b[i]) % MOD;
    }
}

static void mul_x(long long *poly, long long *out) {
    long long overflow = poly[DEGREE - 1];
    for (int i = DEGREE - 1; i > 0; i--) {
        out[i] = poly[i - 1];
    }
    out[0] = 0;
    if (overflow) {
        for (int i = 0; i < DEGREE; i++) {
            out[i] = (out[i] + overflow * FOLD[i]) % MOD;
        }
    }
}

static void poly_mul(long long *a, long long *b, long long *out) {
    long long tmp[DEGREE];
    memset(tmp, 0, sizeof(tmp));
    long long shifted[DEGREE];
    memcpy(shifted, b, sizeof(shifted));
    for (int i = 0; i < DEGREE; i++) {
        if (a[i]) {
            for (int j = 0; j < DEGREE; j++) {
                tmp[j] = (tmp[j] + a[i] * shifted[j]) % MOD;
            }
        }
        long long next_shifted[DEGREE];
        mul_x(shifted, next_shifted);
        memcpy(shifted, next_shifted, sizeof(shifted));
    }
    memcpy(out, tmp, sizeof(tmp));
}

static void poly_pow(long long *base, long long exponent, long long *result) {
    memcpy(result, ONE, sizeof(ONE));
    long long power[DEGREE];
    memcpy(power, base, sizeof(power));
    while (exponent) {
        if (exponent & 1) {
            long long tmp[DEGREE];
            poly_mul(power, result, tmp);
            memcpy(result, tmp, sizeof(tmp));
        }
        exponent >>= 1;
        if (exponent) {
            long long tmp[DEGREE];
            poly_mul(power, power, tmp);
            memcpy(power, tmp, sizeof(tmp));
        }
    }
}

static void poly_compose(long long *outer, long long *inner, long long *result) {
    memcpy(result, ZERO, sizeof(ZERO));
    for (int i = DEGREE - 1; i >= 0; i--) {
        long long tmp[DEGREE];
        poly_mul(result, inner, tmp);
        memcpy(result, tmp, sizeof(tmp));
        result[0] = (result[0] + outer[i]) % MOD;
    }
}

static void iterate_composition(long long *poly, long long count, long long *result) {
    memcpy(result, X, sizeof(X));
    long long power[DEGREE];
    memcpy(power, poly, sizeof(power));
    while (count) {
        if (count & 1) {
            long long tmp[DEGREE];
            poly_compose(power, result, tmp);
            memcpy(result, tmp, sizeof(tmp));
        }
        count >>= 1;
        if (count) {
            long long tmp[DEGREE];
            poly_compose(power, power, tmp);
            memcpy(power, tmp, sizeof(tmp));
        }
    }
}

static void d1(long long length, long long *out) {
    long long base[DEGREE];
    poly_pow(X, length, base);
    long long base_x[DEGREE];
    mul_x(base, base_x);
    poly_add(base, base_x, out);
}

static void d2(long long count, long long *poly, long long *out) {
    long long iter[DEGREE];
    iterate_composition(poly, count, iter);
    long long poly_x[DEGREE];
    mul_x(poly, poly_x);
    poly_compose(iter, poly_x, out);
}

static long long evaluate(long long *poly, long long x) {
    long long result = 0;
    for (int i = DEGREE - 1; i >= 0; i--) {
        result = (result * x + poly[i]) % MOD;
    }
    return result;
}

long long p910_native(void) {
    init_globals();

    long long nesting = 12, middle_count = 345678, length = 9012345;
    long long x_val = 678, addend = 90;

    long long base[DEGREE];
    d1(length, base);

    long long d2val[DEGREE];
    d2(middle_count, base, d2val);

    long long current[DEGREE];
    poly_compose(base, d2val, current);

    for (long long i = 0; i < nesting; i++) {
        long long tmp[DEGREE];
        d2(middle_count, current, tmp);
        memcpy(current, tmp, sizeof(tmp));
    }

    long long val = evaluate(current, x_val);
    return (val + addend) % MOD;
}
