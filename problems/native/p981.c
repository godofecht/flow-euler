// Project Euler 981
// Port of the Python reference solver to C.
// Counts equivalence classes of words over {x,y,z} with cube counts
// (X,Y,Z) where X,Y,Z are cubes below 88^3, weighted by the
// canonical-product sign and q-multinomial at q=-1, summed mod 888888883.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define MOD 888888883LL

static long long mod_pow(long long base, long long exp, long long m) {
    long long result = 1 % m;
    base %= m;
    if (base < 0) base += m;
    while (exp > 0) {
        if (exp & 1) result = result * base % m;
        base = base * base % m;
        exp >>= 1;
    }
    return result;
}

long long p981_native(void) {
    int NC = 88; // cubes for i in range(88): 0..87
    long long *cubes = (long long*)malloc(NC * sizeof(long long));
    for (int i = 0; i < NC; i++) {
        long long ii = (long long)i;
        cubes[i] = ii * ii * ii;
    }

    long long max_n = 3 * cubes[NC - 1]; // 3 * 87^3 = 1975509

    // factorials and inverse factorials mod MOD
    long long *fact = (long long*)malloc((max_n + 1) * sizeof(long long));
    long long *invfact = (long long*)malloc((max_n + 1) * sizeof(long long));

    fact[0] = 1 % MOD;
    for (long long i = 1; i <= max_n; i++) {
        fact[i] = fact[i - 1] * (i % MOD) % MOD;
    }
    invfact[max_n] = mod_pow(fact[max_n], MOD - 2, MOD);
    for (long long i = max_n; i >= 1; i--) {
        invfact[i - 1] = invfact[i] * (i % MOD) % MOD;
    }

    long long inv2 = (MOD + 1) / 2; // inverse of 2 mod odd MOD

    // per-cube helpers
    long long *halves = (long long*)malloc(NC * sizeof(long long));
    long long *invf = (long long*)malloc(NC * sizeof(long long));
    int *par = (int*)malloc(NC * sizeof(int));
    for (int i = 0; i < NC; i++) {
        halves[i] = cubes[i] >> 1;
        invf[i] = invfact[cubes[i]];
        par[i] = i & 1; // same parity as cube (cube parity == i parity)
    }

    long long total_sum = 0;

    for (int ai = 0; ai < NC; ai++) {
        long long X = cubes[ai];
        long long hx = halves[ai];
        long long invX = invf[ai];
        int px = par[ai];
        for (int bj = 0; bj < NC; bj++) {
            long long Y = cubes[bj];
            long long hy = halves[bj];
            long long invY = invf[bj];
            int py = par[bj];
            if (px != py) continue;
            for (int ck = 0; ck < NC; ck++) {
                if (py != par[ck]) continue;

                long long Z = cubes[ck];
                long long hz = halves[ck];
                long long invZ = invf[ck];

                long long n = X + Y + Z;

                // T = multinomial(n; X,Y,Z) mod MOD
                long long T = fact[n];
                T = T * invX % MOD;
                T = T * invY % MOD;
                T = T * invZ % MOD;

                // D = (E - O) = q-multinomial at q=-1
                long long D;
                if ((n & 1) == 0 && (X & 1) == 1) {
                    D = 0;
                } else {
                    // D1 = comb_mod(n>>1, hx)
                    long long d1n = n >> 1;
                    long long D1;
                    if (hx < 0 || hx > d1n) {
                        D1 = 0;
                    } else {
                        D1 = fact[d1n] * invfact[hx] % MOD * invfact[d1n - hx] % MOD;
                    }
                    long long n2 = n - X; // Y + Z
                    if ((n2 & 1) == 0 && (Y & 1) == 1) {
                        D = 0;
                    } else {
                        long long d2n = n2 >> 1;
                        long long D2;
                        if (hy < 0 || hy > d2n) {
                            D2 = 0;
                        } else {
                            D2 = fact[d2n] * invfact[hy] % MOD * invfact[d2n - hy] % MOD;
                        }
                        D = D1 * D2 % MOD;
                    }
                }

                // sign: +1 if (hx+hy+hz) even else -1
                long long Nmod;
                if (((hx + hy + hz) & 1) == 0) {
                    Nmod = (T + D) % MOD;
                } else {
                    Nmod = (T - D) % MOD;
                    if (Nmod < 0) Nmod += MOD;
                }

                Nmod = Nmod * inv2 % MOD;
                total_sum += Nmod;
                if (total_sum >= MOD) total_sum -= MOD;
            }
        }
    }

    long long answer = total_sum % MOD;

    free(cubes);
    free(fact);
    free(invfact);
    free(halves);
    free(invf);
    free(par);

    return answer;
}
