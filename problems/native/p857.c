/* Project Euler 857: Beautiful Graphs.
   Computes G(10^7) mod (10^9+7).

   a[s] = number of 2-colourings of the edges of K_s with no
   monochromatic triangle (brute force for s <= 5).

   G(n) satisfies the recurrence (ordinary generating function):
       f_n = sum_{s=1}^{5} (a_s / s!) * f_{n-s},  f_0 = 1
   and G(n) = n! * f_n  (mod MOD).
*/

#include <stdint.h>
#include <stdio.h>

typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

#define MOD 1000000007LL
#define TARGET_N 10000000LL

static i64 mod_pow(i64 a, i64 e, i64 mod) {
    i64 r = 1 % mod;
    a %= mod;
    if (a < 0) a += mod;
    while (e > 0) {
        if (e & 1) r = (i64)((i128)r * a % mod);
        a = (i64)((i128)a * a % mod);
        e >>= 1;
    }
    return r;
}

static i64 mod_inv(i64 a, i64 mod) { return mod_pow(a, mod - 2, mod); }

/* Count 2-colourings of edges of K_k with no monochromatic triangle. */
static i64 count_no_mono_triangle(int k) {
    if (k <= 1) return 1;

    /* edge index: (i,j) -> idx, i < j */
    int edge[6][6];
    int m = 0;
    for (int i = 0; i < k; i++)
        for (int j = 0; j < k; j++)
            edge[i][j] = -1;
    for (int i = 0; i < k; i++) {
        for (int j = i + 1; j < k; j++) {
            edge[i][j] = m;
            m++;
        }
    }

    /* collect triangles as triples of edge indices */
    int tris[20][3];
    int nt = 0;
    for (int i = 0; i < k; i++) {
        for (int j = i + 1; j < k; j++) {
            for (int l = j + 1; l < k; l++) {
                tris[nt][0] = edge[i][j];
                tris[nt][1] = edge[i][l];
                tris[nt][2] = edge[j][l];
                nt++;
            }
        }
    }

    i64 good = 0;
    u64 total = (u64)1 << m;
    for (u64 mask = 0; mask < total; mask++) {
        int ok = 1;
        for (int t = 0; t < nt && ok; t++) {
            int x = (int)((mask >> tris[t][0]) & 1);
            int y = (int)((mask >> tris[t][1]) & 1);
            int z = (int)((mask >> tris[t][2]) & 1);
            if (x == y && x == z) ok = 0;
        }
        if (ok) good++;
    }
    return good;
}

static i64 G_mod(i64 n) {
    if (n == 0) return 1;

    i64 a[6];
    a[0] = 0;
    for (int s = 1; s <= 5; s++) a[s] = count_no_mono_triangle(s);

    /* coeff[s] = a[s] / s!  (mod MOD) */
    i64 coeff[6];
    i64 fact_small = 1;
    for (int s = 1; s <= 5; s++) {
        fact_small = (i64)((i128)fact_small * s % MOD);
        coeff[s] = (i64)((i128)a[s] * mod_inv(fact_small, MOD) % MOD);
    }

    i64 c1 = coeff[1], c2 = coeff[2], c3 = coeff[3], c4 = coeff[4], c5 = coeff[5];

    /* rolling recurrence: f1 = f_0 = 1, f2..f5 = 0 (lagged) */
    i64 f1 = 1, f2 = 0, f3 = 0, f4 = 0, f5 = 0;

    i64 fact = 1; /* n! mod MOD */
    for (i64 i = 1; i <= n; i++) {
        i64 fn = (i64)((i128)c1 * f1 % MOD
                     + (i128)c2 * f2 % MOD
                     + (i128)c3 * f3 % MOD
                     + (i128)c4 * f4 % MOD
                     + (i128)c5 * f5 % MOD) % MOD;
        f5 = f4; f4 = f3; f3 = f2; f2 = f1; f1 = fn;
        fact = (i64)((i128)fact * i % MOD);
    }

    return (i64)((i128)fact * f1 % MOD);
}

long long p857_native(void) {
    return G_mod(TARGET_N);
}
