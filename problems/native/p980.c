// Project Euler 980
// Count ordered pairs (i, j) with 0 <= i, j < N such that the
// concatenation c(i)c(j) of 50-letter {x,y,z} strings is neutral.
// Neutrality maps to the quaternion group Q8 product being +1 under
// x -> i, y -> j, z -> -k.
// Port of the Python reference solver to C.
#include <stdint.h>
#include <stdio.h>

static void build_mul_table_q8(int mul[8][8]) {
    /* Elements 0..7: 0:1, 1:i, 2:j, 3:k, 4:-1, 5:-i, 6:-j, 7:-k */
    int base[4][4];
    int sgn[4][4];

    for (int a = 0; a < 4; a++) {
        for (int b = 0; b < 4; b++) {
            if (a == 0) {
                base[a][b] = b;
                sgn[a][b] = 1;
            } else if (b == 0) {
                base[a][b] = a;
                sgn[a][b] = 1;
            } else if (a == b) {
                base[a][b] = 0;
                sgn[a][b] = -1;
            } else {
                if (a == 1 && b == 2) { base[a][b] = 3; sgn[a][b] = 1; }
                else if (a == 2 && b == 3) { base[a][b] = 1; sgn[a][b] = 1; }
                else if (a == 3 && b == 1) { base[a][b] = 2; sgn[a][b] = 1; }
                else if (a == 2 && b == 1) { base[a][b] = 3; sgn[a][b] = -1; }
                else if (a == 3 && b == 2) { base[a][b] = 1; sgn[a][b] = -1; }
                else if (a == 1 && b == 3) { base[a][b] = 2; sgn[a][b] = -1; }
            }
        }
    }

    for (int A = 0; A < 8; A++) {
        int sa = (A < 4) ? 1 : -1;
        int a = A & 3;
        for (int B = 0; B < 8; B++) {
            int sb = (B < 4) ? 1 : -1;
            int b = B & 3;
            int s = sa * sb * sgn[a][b];
            int c = base[a][b];
            mul[A][B] = (s == 1) ? c : (c ^ 4);
        }
    }
}

long long p980_native(void) {
    int mul[8][8];
    build_mul_table_q8(mul);

    /* x -> i (1), y -> j (2), z -> -k (7) */
    int gen[3] = {1, 2, 7};

    /* Right-multiply table: R[v*3 + b] = mul[v][gen[b]] */
    int R[8 * 3];
    for (int v = 0; v < 8; v++)
        for (int b = 0; b < 3; b++)
            R[v * 3 + b] = mul[v][gen[b]];

    /* Inverses in Q8 */
    int inv[8];
    for (int e = 0; e < 8; e++) {
        for (int f = 0; f < 8; f++) {
            if (mul[e][f] == 0 && mul[f][e] == 0) {
                inv[e] = f;
                break;
            }
        }
    }

    const long long MOD = 888888883LL;
    const long long MULT = 8888LL;
    long long a = 88888888LL; /* a_0 */

    long long N = 1000000LL;
    long long counts[8];
    for (int i = 0; i < 8; i++) counts[i] = 0;

    for (long long i = 0; i < N; i++) {
        int v = 0; /* quaternion product for this 50-letter block */
        for (int j = 0; j < 50; j++) {
            v = R[v * 3 + (int)(a % 3)];
            a = (a * MULT) % MOD;
        }
        counts[v]++;
    }

    long long total = 0;
    for (int e = 0; e < 8; e++) {
        total += counts[e] * counts[inv[e]];
    }
    return total;
}
