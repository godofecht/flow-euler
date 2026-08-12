// Project Euler 895: Gold & Silver Coin Game II.
// Count ordered triples of stacks that are fair and balanced.
// G(9898) mod 989898989.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

static i64 mulmod(i64 a, i64 b, i64 m) {
    return (i64)((i128)a * b % m);
}

static i64 powmod(i64 a, i64 e, i64 m) {
    i64 r = 1 % m; a %= m; if (a < 0) a += m;
    while (e > 0) {
        if (e & 1) r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}

static i64 modinv(i64 a, i64 m) {
    a %= m; if (a < 0) a += m;
    i64 x0 = 1, x1 = 0, aa = a, mm = m;
    while (mm) {
        i64 q = aa / mm;
        i64 t = aa - q * mm; aa = mm; mm = t;
        t = x0 - q * x1; x0 = x1; x1 = t;
    }
    if (aa != 1) return -1;
    x0 %= m; if (x0 < 0) x0 += m;
    return x0;
}

static i64 ceil_div(i64 n, i64 d) {
    return -((-n) / d);
}

// C(S+2, 2) for S >= 0, else 0
static i64 tri(i64 S) {
    if (S < 0) return 0;
    i64 t = S + 2;
    return t * (t - 1) / 2;
}

// Count integer solutions y_i in [0, n_i-1] with sum = target.
// Inclusion-exclusion over which variables exceed upper bound.
static i64 bounded_sum_count(i64 n0, i64 n1, i64 n2, i64 target) {
    i64 total = 0;
    for (int mask = 0; mask < 8; mask++) {
        i64 sub = 0;
        int bits = 0;
        if (mask & 1) { sub += n0; bits++; }
        if (mask & 2) { sub += n1; bits++; }
        if (mask & 4) { sub += n2; bits++; }
        i64 val = tri(target - sub);
        if (bits & 1) total -= val;
        else total += val;
    }
    return total;
}

// ---- Fast modular solver ----

static i64 MOD;
static i64 inv2;
static i64 *pow2_arr, *invpow2_arr;
static i64 *P0, *P1, *P2;

static void interval_sums(int l, int r, i64 *s0, i64 *s1, i64 *s2) {
    if (l > r) { *s0 = 0; *s1 = 0; *s2 = 0; return; }
    *s0 = (P0[r] - P0[l - 1]) % MOD; if (*s0 < 0) *s0 += MOD;
    *s1 = (P1[r] - P1[l - 1]) % MOD; if (*s1 < 0) *s1 += MOD;
    *s2 = (P2[r] - P2[l - 1]) % MOD; if (*s2 < 0) *s2 += MOD;
}

static i64 sum_F_linear(i64 alpha, i64 beta, int l, int r) {
    if (l > r) return 0;
    i64 s0, s1, s2;
    interval_sums(l, r, &s0, &s1, &s2);
    i64 a_mod = alpha % MOD; if (a_mod < 0) a_mod += MOD;
    i64 b_mod = beta % MOD; if (b_mod < 0) b_mod += MOD;

    i64 term2 = mulmod(a_mod, a_mod, MOD);
    i64 term1 = mulmod(a_mod, (2 * b_mod + 3) % MOD, MOD);
    i64 term0 = (b_mod * b_mod + 3 * b_mod + 2) % MOD;

    i64 res = (mulmod(term2, s2, MOD) + mulmod(term1, s1, MOD) + mulmod(term0, s0, MOD)) % MOD;
    return mulmod(res, inv2, MOD);
}

static const i64 C2[3] = {1, 2, 1};

static i64 G_pq(int b, int s, i64 p, i64 q) {
    int Amax = b - 1;
    i64 total = 0;
    for (int ca = 0; ca <= 2; ca++) {
        i64 mult = C2[ca];
        for (int cb = 0; cb <= 1; cb++) {
            i64 sign = ((ca + cb) & 1) ? -1 : 1;
            i64 coeff = sign * mult;
            i64 alpha = p - ca;
            i64 beta = (q - cb) * b - s;

            int l, r;
            if (alpha == 0) {
                if (beta < 0) continue;
                l = 1; r = Amax;
            } else if (alpha > 0) {
                l = (int)ceil_div(-beta, alpha);
                if (l < 1) l = 1;
                r = Amax;
                if (l > r) continue;
            } else {
                r = (int)(beta / (-alpha));
                if (r > Amax) r = Amax;
                l = 1;
                if (r < l) continue;
            }

            i64 val = sum_F_linear(alpha, beta, l, r);
            total = (total + coeff * val) % MOD;
        }
    }
    return total % MOD;
}

static void base_weighted(int b, int s, i64 *base) {
    i64 G[3][2];
    for (int p = 0; p <= 2; p++)
        for (int q = 0; q <= 1; q++)
            G[p][q] = G_pq(b, s, p, q);

    for (int r = 0; r <= 3; r++) {
        i64 acc = 0;
        for (int nb = 0; nb <= 1; nb++) {
            int ra = r - nb;
            if (ra >= 0 && ra <= 2) {
                i64 mult_sign = 3 * C2[ra];
                acc = (acc + mulmod(mult_sign, G[ra][nb], MOD)) % MOD;
            }
        }
        base[r] = mulmod(acc, pow2_arr[b - 1], MOD);
    }
}

static i64 G_mod(int m, i64 mod) {
    MOD = mod;
    inv2 = modinv(2, MOD);

    pow2_arr = (i64 *)malloc((m + 1) * sizeof(i64));
    invpow2_arr = (i64 *)malloc((m + 1) * sizeof(i64));
    P0 = (i64 *)malloc((m + 1) * sizeof(i64));
    P1 = (i64 *)malloc((m + 1) * sizeof(i64));
    P2 = (i64 *)malloc((m + 1) * sizeof(i64));

    pow2_arr[0] = 1;
    for (int i = 1; i <= m; i++)
        pow2_arr[i] = mulmod(pow2_arr[i - 1], 2, MOD);

    invpow2_arr[0] = 1;
    invpow2_arr[1] = inv2 % MOD;
    for (int i = 2; i <= m; i++)
        invpow2_arr[i] = mulmod(invpow2_arr[i - 1], inv2, MOD);

    P0[0] = 0; P1[0] = 0; P2[0] = 0;
    for (int a = 1; a <= m; a++) {
        i64 w = invpow2_arr[a];
        P0[a] = (P0[a - 1] + w) % MOD;
        P1[a] = (P1[a - 1] + mulmod(a, w, MOD)) % MOD;
        P2[a] = (P2[a - 1] + mulmod(mulmod(a, a, MOD), w, MOD)) % MOD;
    }

    // Case 0: 3 monochrome
    i64 case0 = mulmod(3 * m, (m - 1), MOD);

    // Case 2: two mixed + one monochrome
    i64 case2 = 0;
    for (int t = 1; t < m; t++) {
        int n = m - t;
        i64 term = mulmod(pow2_arr[t - 1], mulmod(n, (n - 1), MOD), MOD);
        case2 = (case2 + term) % MOD;
    }
    case2 = mulmod(case2, 6, MOD);

    // Case 3: three mixed
    i64 case3 = 0;

    // Carry DP for u-bit core
    i64 *cur0 = (i64 *)malloc((m + 1) * sizeof(i64));
    i64 *cur1 = (i64 *)malloc((m + 1) * sizeof(i64));
    i64 *nxt0 = (i64 *)malloc((m + 1) * sizeof(i64));
    i64 *nxt1 = (i64 *)malloc((m + 1) * sizeof(i64));

    cur0[0] = 1;  // u=1
    cur1[0] = 1;

    for (int u = 1; u <= m - 1; u++) {
        if (u <= m - 2) {
            int b = m - u;

            i64 base_s1[4], base_s2[4];
            base_weighted(b, 1, base_s1);
            base_weighted(b, 2, base_s2);

            // s=1, f=0 and s=2, f=1
            for (int si = 0; si < 2; si++) {
                int s = si + 1;
                int f = si;
                i64 *base = (si == 0) ? base_s1 : base_s2;
                for (int r = 1; r <= 3; r++) {
                    i64 Wtarget = s - r;
                    i64 num = Wtarget + u + 1 - 4 * f;
                    if (num & 1) continue;
                    i64 C = num / 2;
                    if (C < 0 || C > u - 1) continue;
                    i64 numerator_high = (f == 0) ? cur0[C] : cur1[C];
                    case3 = (case3 + mulmod(numerator_high, base[r], MOD)) % MOD;
                }
            }
        }

        // Update DP to u+1
        memset(nxt0, 0, (u + 1) * sizeof(i64));
        memset(nxt1, 0, (u + 1) * sizeof(i64));
        nxt0[0] = mulmod(3, cur0[0], MOD);
        nxt1[0] = cur0[0] % MOD;
        for (int c = 1; c < u; c++) {
            nxt0[c] = (mulmod(3, cur0[c], MOD) + cur1[c - 1]) % MOD;
            nxt1[c] = (cur0[c] + mulmod(3, cur1[c - 1], MOD)) % MOD;
        }
        nxt0[u] = cur1[u - 1] % MOD;
        nxt1[u] = mulmod(3, cur1[u - 1], MOD);

        memcpy(cur0, nxt0, (u + 1) * sizeof(i64));
        memcpy(cur1, nxt1, (u + 1) * sizeof(i64));
    }

    free(cur0); free(cur1); free(nxt0); free(nxt1);
    free(pow2_arr); free(invpow2_arr);
    free(P0); free(P1); free(P2);

    return (case0 + case2 + case3) % MOD;
}

long long p895_native(void) {
    int m = 9898;
    i64 mod = 989898989;
    return (long long)G_mod(m, mod);
}
