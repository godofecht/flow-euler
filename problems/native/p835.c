#include <stdio.h>
#include <math.h>

/* Project Euler 835 - Supernatural Triangles.
 *
 * S(N) = sum of perimeters of all distinct supernatural triangles with P <= N.
 * Compute S(10^(10^10)) mod 1234567891.
 *
 * Family B: c = b+1, a odd. P = t(t+1), odd t >= 3, t <= 10^M - 1, M = E/2.
 * Family A: legs consecutive. Pell x^2 - 2y^2 = -1. P_k recurrence.
 * Overlap: (3,4,5) perimeter 12 in both families.
 */

typedef long long i64;
typedef __int128 i128;

static const i64 MOD = 1234567891LL;

static i64 mmul(i64 a, i64 b) {
    return (i64)((i128)a * b % MOD);
}

static i64 mpow(i64 base, i64 exp) {
    i64 result = 1 % MOD;
    base %= MOD;
    if (base < 0) base += MOD;
    while (exp > 0) {
        if (exp & 1) result = mmul(result, base);
        base = mmul(base, base);
        exp >>= 1;
    }
    return result;
}

static i64 egcd(i64 a, i64 b, i64 *x, i64 *y) {
    if (b == 0) { *x = 1; *y = 0; return a; }
    i64 x1, y1;
    i64 g = egcd(b, a % b, &x1, &y1);
    *x = y1;
    *y = x1 - (a / b) * y1;
    return g;
}

static i64 modinv(i64 a) {
    i64 x, y;
    a %= MOD;
    if (a < 0) a += MOD;
    egcd(a, MOD, &x, &y);
    return ((x % MOD) + MOD) % MOD;
}

/* ---- Family B ---- */
static i64 sum_family_B(i64 E) {
    /* E is even, E >= 2. M = E/2. */
    i64 M = E / 2;
    /* n = 5 * 10^(M-1) mod MOD */
    i64 n = (5 * mpow(10, M - 1)) % MOD;
    i64 n2 = mmul(n, n);
    i64 n3 = mmul(n2, n);
    i64 INV3 = modinv(3);

    i64 sum_t_all = n2;  /* 1+3+...+(2n-1) = n^2 */
    i64 sum_t2_all = mmul((4 * n3 % MOD - n % MOD + MOD) % MOD, INV3);
    /* Remove t=1 term: 1 + 1 = 2 */
    return (sum_t_all + sum_t2_all - 2 + MOD) % MOD;
}

/* ---- Family A: matrix exponentiation ---- */
static void mat_mul(i64 A[3][3], i64 B[3][3], i64 C[3][3]) {
    i64 tmp[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            i128 s = 0;
            for (int k = 0; k < 3; k++)
                s += (i128)A[i][k] * B[k][j];
            tmp[i][j] = (i64)(s % MOD);
        }
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            C[i][j] = tmp[i][j];
}

static void mat_pow(i64 M[3][3], i64 exp, i64 R[3][3]) {
    i64 result[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    i64 base[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            base[i][j] = M[i][j];
    while (exp > 0) {
        if (exp & 1) mat_mul(result, base, result);
        mat_mul(base, base, base);
        exp >>= 1;
    }
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            R[i][j] = result[i][j];
}

static i64 pell_max_index(i64 E) {
    double alpha = 3.0 + 2.0 * sqrt(2.0);
    double A_coeff = (4.0 + 3.0 * sqrt(2.0)) / 4.0;
    double log10_alpha = log10(alpha);
    double log10_A = log10(A_coeff);

    double x = ((double)E - log10_A) / log10_alpha;
    i64 n = (i64)x;

    /* Adjust for floating point rounding */
    while (n > 0 && (log10_A + (double)n * log10_alpha) >= (double)E) n--;
    while ((log10_A + (double)(n + 1) * log10_alpha) < (double)E) n++;
    return n >= 1 ? n : 0;
}

static i64 sum_family_A(i64 n) {
    if (n <= 0) return 0;
    if (n == 1) return 12 % MOD;

    /* State: [P_k, P_{k-1}, S_k]^T
     * P_{k+1} = 6*P_k - P_{k-1}
     * S_{k+1} = S_k + P_{k+1}
     */
    i64 mat[3][3] = {
        {6, MOD - 1, 0},
        {1, 0, 0},
        {6, MOD - 1, 1}
    };
    i64 P[3][3];
    mat_pow(mat, n - 1, P);

    /* v1 = [12, 2, 12]^T (k=1) */
    i64 v[3] = {12, 2, 12};
    i64 result[3];
    for (int i = 0; i < 3; i++) {
        i128 s = 0;
        for (int j = 0; j < 3; j++)
            s += (i128)P[i][j] * v[j];
        result[i] = (i64)(s % MOD);
    }
    return result[2] % MOD;
}

i64 p835_native(void) {
    i64 E = 10000000000LL; /* 10^10 */

    i64 sumB = sum_family_B(E);
    i64 nA = pell_max_index(E);
    i64 sumA = sum_family_A(nA);

    /* Remove overlap: perimeter 12 in both families */
    i64 ans = (sumA + sumB - 12 + MOD) % MOD;
    return ans;
}
