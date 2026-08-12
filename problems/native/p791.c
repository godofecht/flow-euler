// Project Euler 791: Average and Variance
// Compute S(n) = sum of 2*m over all quadruples (a,b,c,d) with
// 1 <= a <= b <= c <= d <= n whose average equals twice their variance.
// Uses the (U,V,W,sgn) lattice parametrization with prefix-sum aggregation.
// Answer: S(10^8) mod 433494437.
// Port of the Python reference solver to C.
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

static i64 isqrt_i64(i64 n) {
    if (n <= 0) return 0;
    i64 r = (i64)sqrtl((long double)n);
    // correct for any rounding
    while (r > 0 && r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r;
}

long long p791_native(void) {
    i64 n = 100000000LL;  // 10^8
    i64 mod = 433494437LL;
    i64 N2 = 2 * n;
    i64 R = isqrt_i64(2 * n) + 2;

    // Prefix sums
    i64 *s2 = calloc(R + 1, sizeof(i64));
    i64 *ps2 = calloc(R + 1, sizeof(i64));
    i64 *p1 = calloc(R + 1, sizeof(i64));
    i64 *p2 = calloc(R + 1, sizeof(i64));
    i64 *p3 = calloc(R + 1, sizeof(i64));
    if (!s2 || !ps2 || !p1 || !p2 || !p3) return 0;

    for (i64 i = 1; i <= R; i++) {
        i64 ii = i * i;
        s2[i] = s2[i - 1] + ii;
        ps2[i] = ps2[i - 1] + s2[i];
        p1[i] = p1[i - 1] + i;
        p2[i] = p2[i - 1] + ii;
        p3[i] = p3[i - 1] + ii * i;
    }

    i64 total = 0;
    const i64 LIM = (i64)1 << 62;  // use 2^62 to stay safe in signed range

    // Small U cases (U=0,1): brute force
    for (int U = 0; U < 2; U++) {
        for (int V = 0; V <= U; V++) {
            for (int W = 0; W <= V; W++) {
                for (int sgn = -1; sgn <= 1; sgn += 2) {
                    if (W == 0 && sgn == -1) continue;
                    i64 u = -U;
                    i64 v = -V;
                    i64 w = (i64)sgn * W;
                    i64 m = (i64)U * U + (i64)V * V + (i64)W * W;
                    if (m == 0) continue;
                    i64 a = (m + u + v + w) / 2;
                    i64 b = (m + u - v - w) / 2;
                    i64 c = (m - u + v - w) / 2;
                    i64 d = (m - u - v + w) / 2;
                    if (1 <= a && a <= b && b <= c && c <= d && d <= n) {
                        i64 v_add = 2 * m;
                        total += v_add % mod;
                        if (total >= LIM) total %= mod;
                    }
                }
            }
        }
    }

    // Main enumeration
    for (i64 U = 2; U <= R; U++) {
        i64 U2 = U * U;
        i64 rem = N2 - U2 - U;
        if (rem < 0) break;

        i64 Wmax0 = isqrt_i64(rem / 2);
        if (Wmax0 > U) Wmax0 = U;

        i64 T = N2 - 2 * U2 - 2 * U;

        for (int sgn = -1; sgn <= 1; sgn += 2) {
            i64 startW = (sgn == 1) ? 0 : 1;
            if (startW > Wmax0) continue;

            i64 Wfull;
            if (T < 0) {
                Wfull = -1;
            } else {
                i64 rt = isqrt_i64(1 + 4 * T);
                if (sgn == 1) {
                    Wfull = (rt - 1) / 2;
                } else {
                    Wfull = (rt + 1) / 2;
                }
                if (Wfull > Wmax0) Wfull = Wmax0;
            }

            // Full segment: W in [startW..Wfull], V in [W..U]
            if (Wfull >= startW) {
                i64 A = startW;
                i64 B = Wfull;
                i64 num = B - A + 1;

                i64 sumW = p1[B] - (A ? p1[A - 1] : 0);
                i64 sumW2 = p2[B] - (A ? p2[A - 1] : 0);
                i64 sumW3 = p3[B] - (A ? p3[A - 1] : 0);

                i64 sumCnt = num * (U + 1) - sumW;
                i64 sumCntW2 = (U + 1) * sumW2 - sumW3;

                i64 sumPrefix;
                if (B == 0) {
                    sumPrefix = 0;
                } else {
                    i64 lo = A - 1;
                    i64 hi = B - 1;
                    if (lo < 0) lo = 0;
                    sumPrefix = ps2[hi] - (lo > 0 ? ps2[lo - 1] : 0);
                }

                i64 sumSumV2 = num * s2[U] - sumPrefix;

                // contrib = 2*U2*sumCnt + 2*sumCntW2 + 2*sumSumV2
                // Use i128 to avoid overflow, then reduce mod
                i128 contrib = (i128)2 * U2 * sumCnt + (i128)2 * sumCntW2 + (i128)2 * sumSumV2;
                i64 cmod = (i64)(contrib % mod);
                if (cmod < 0) cmod += mod;
                total += cmod;
                if (total >= LIM) total %= mod;
            }

            // Tail segment
            i64 tail_start = Wfull + 1;
            if (tail_start < startW) tail_start = startW;
            if (tail_start > Wmax0) continue;

            i64 base_D = 4 * N2 + 1 - 4 * (U2 + U);

            for (i64 W = tail_start; W <= Wmax0; W++) {
                i64 D = base_D - 4 * (W * W + (i64)sgn * W);
                if (D < 0) break;

                i64 Vmax = (isqrt_i64(D) - 1) / 2;
                if (Vmax > U) Vmax = U;
                if (Vmax < W) continue;

                i64 cnt = Vmax - W + 1;
                i64 sumV2 = s2[Vmax] - (W ? s2[W - 1] : 0);

                i64 W2 = W * W;
                i128 contrib = (i128)cnt * 2 * (U2 + W2) + (i128)2 * sumV2;
                i64 cmod = (i64)(contrib % mod);
                if (cmod < 0) cmod += mod;
                total += cmod;
                if (total >= LIM) total %= mod;
            }
        }
    }

    free(s2); free(ps2); free(p1); free(p2); free(p3);
    return total % mod;
}
