// Project Euler 803: Pseudorandom Sequence
// Find the minimum index n such that the rand48-based pseudorandom sequence
// starting from the unique seed for "PuzzleOne" produces "LuckyText".
// Uses 48-bit LCG state splitting, CRT-based seed recovery, and
// 2-adic discrete logarithm.
// Port of the Python reference solver to C.
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef uint64_t u64;
typedef __uint128_t u128;

#define A_VAL 25214903917ULL
#define C_VAL 11ULL
#define MOD48 (1ULL << 48)
#define MASK48 (MOD48 - 1)
#define MOD24 (1ULL << 24)
#define MASK24 (MOD24 - 1)
#define MOD18 (1ULL << 18)
#define MASK18 (MOD18 - 1)
#define INV9_MOD13 3
#define ORDER_EXP 46

// ---- mod-2^48 arithmetic ----

static u64 mul48(u64 a, u64 b) {
    return (u64)((u128)a * b & MASK48);
}

static u64 powmod48(u64 base, u64 exp) {
    u64 result = 1;
    base &= MASK48;
    while (exp > 0) {
        if (exp & 1) result = mul48(result, base);
        base = mul48(base, base);
        exp >>= 1;
    }
    return result;
}

static u64 invmod48(u64 x) {
    x &= MASK48;
    // x is odd: x^2 ≡ 1 mod 8, so inv=x is correct mod 8.
    // Newton: inv = inv*(2 - x*inv) mod 2^48, doubling precision each step.
    u64 inv = x;
    for (int i = 0; i < 5; i++) {  // 3->6->12->24->48->96 (extra is harmless)
        u64 t = mul48(x, inv);
        u64 two_minus = (0 - t + 2) & MASK48;  // (2 - t) mod 2^48
        inv = mul48(inv, two_minus);
    }
    return inv;
}

static u64 step48(u64 x) {
    return (mul48(A_VAL, x) + C_VAL) & MASK48;
}

// ---- character mapping ----

static int ch_to_val(char ch) {
    int o = (unsigned char)ch;
    if (o >= 97 && o <= 122) return o - 97;       // a-z -> 0-25
    return o - 65 + 26;                            // A-Z -> 26-51
}

static char val_to_ch(int v) {
    if (v < 26) return (char)(97 + v);
    return (char)(65 + (v - 26));
}

static int b_from_a(u64 a) {
    return (int)((a >> 16) % 52);
}

// Check if seed a0 produces the given prefix of values
static int check_prefix(u64 a0, const int *vals, int L) {
    u64 a = a0 & MASK48;
    for (int i = 0; i < L; i++) {
        if (b_from_a(a) != vals[i]) return 0;
        a = step48(a);
    }
    return 1;
}

// ---- u0 candidates: brute-force u0 in [0,2^18) satisfying mod-4 constraints ----

static int u0_candidates(const int *pattern_vals, int L, u64 *out, int max_out) {
    int count = 0;
    for (u64 u0 = 0; u0 < MOD18; u0++) {
        u64 u = u0;
        int ok = 1;
        for (int i = 0; i < L; i++) {
            if (((u >> 16) & 3) != (pattern_vals[i] & 3)) {
                ok = 0;
                break;
            }
            u = (mul48(A_VAL, u) + C_VAL) & MASK18;
        }
        if (ok) {
            if (count < max_out) out[count] = u0;
            count++;
        }
    }
    return count;
}

// ---- solve y0 for residues ----
// Given y_{n+1} = (A*y_n + carries24[n]) mod 2^24, and y_n mod 13 = residues13[n],
// find all y0 in [0, 2^24) satisfying the constraints.

static int solve_y0(u64 *carries24, int *residues13, int L,
                    u64 *out, int max_out) {
    if (L == 0) return 0;
    int r0 = residues13[0];
    int count = 0;

    if (L == 1) {
        for (u64 y0 = r0; y0 < MOD24; y0 += 13) {
            if (count < max_out) out[count] = y0;
            count++;
        }
        return count;
    }

    int r1 = residues13[1];
    u64 a24 = A_VAL & MASK24;

    if (L == 2) {
        u64 y0 = r0;
        u64 y1 = (mul48(A_VAL, y0) + carries24[0]) & MASK24;
        u64 delta1 = (13 * a24) & MASK24;
        for (y0 = r0; y0 < MOD24; y0 += 13) {
            if ((y1 % 13) == (u64)r1) {
                if (count < max_out) out[count] = y0;
                count++;
            }
            y1 = (y1 + delta1) & MASK24;
        }
        return count;
    }

    // L >= 3
    int r2 = residues13[2];
    u64 y0 = r0;
    u64 y1 = (mul48(A_VAL, y0) + carries24[0]) & MASK24;
    u64 y2 = (mul48(A_VAL, y1) + carries24[1]) & MASK24;
    u64 delta1 = (13 * a24) & MASK24;
    u64 delta2 = mul48(delta1, a24);

    for (y0 = r0; y0 < MOD24; y0 += 13) {
        if ((y1 % 13) == (u64)r1 && (y2 % 13) == (u64)r2) {
            // full verify from y2 onward
            u64 y = y2;
            int ok = 1;
            for (int i = 2; i < L - 1; i++) {
                y = (mul48(A_VAL, y) + carries24[i]) & MASK24;
                if ((y % 13) != (u64)residues13[i + 1]) {
                    ok = 0;
                    break;
                }
            }
            if (ok) {
                if (count < max_out) out[count] = y0;
                count++;
            }
        }
        y1 = (y1 + delta1) & MASK24;
        y2 = (y2 + delta2) & MASK24;
    }
    return count;
}

// ---- solve states for pattern ----

static int solve_states(const char *pattern, u64 *states_out, int max_states) {
    int vals[64];
    int L = (int)strlen(pattern);
    for (int i = 0; i < L; i++) vals[i] = ch_to_val(pattern[i]);

    u64 us[1024];
    int n_us = u0_candidates(vals, L, us, 1024);

    int n_states = 0;

    u64 u_list[64], k_list[64];
    u64 carries24[64];
    int t_list[64], residues13[64];
    u64 y0s[4096];

    for (int ui = 0; ui < n_us; ui++) {
        u64 u0 = us[ui];
        // Precompute u_list and k_list
        u64 u = u0;
        for (int i = 0; i < L; i++) {
            u_list[i] = u;
            u128 nxt = (u128)A_VAL * u + C_VAL;
            if (i < L - 1) k_list[i] = (u64)(nxt >> 18);
            u = (u64)(nxt & MASK18);
        }

        // Enumerate w0 (6 bits)
        for (int w0 = 0; w0 < 64; w0++) {
            u64 w = w0;
            for (int i = 0; i < L; i++) {
                t_list[i] = (int)(((u_list[i] >> 16) & 3) + (w << 2));
                if (i < L - 1) {
                    carries24[i] = (k_list[i] + A_VAL * w) >> 6;
                    w = (A_VAL * w + k_list[i]) & 63;
                }
            }

            // Compute residues13
            for (int i = 0; i < L; i++) {
                int t_mod13 = t_list[i] % 13;
                if (t_mod13 < 0) t_mod13 += 13;
                int diff = (vals[i] - t_mod13) % 13;
                if (diff < 0) diff += 13;
                residues13[i] = (INV9_MOD13 * diff) % 13;
            }

            // Solve y0
            int n_y0 = solve_y0(carries24, residues13, L, y0s, 4096);
            for (int yi = 0; yi < n_y0; yi++) {
                u64 y0 = y0s[yi];
                u64 x0 = u0 + ((u64)w0 << 18);
                u64 a0 = x0 + (y0 << 24);
                if (check_prefix(a0, vals, L)) {
                    if (n_states < max_states) states_out[n_states] = a0;
                    n_states++;
                }
            }
        }
    }
    return n_states;
}

// ---- 2-adic discrete log ----

static u64 G_ORDER2;
static u64 INV_POWS_2I[ORDER_EXP];

static void precompute_dlog(void) {
    G_ORDER2 = powmod48(A_VAL, 1ULL << (ORDER_EXP - 1));
    for (int i = 0; i < ORDER_EXP; i++) {
        u64 p = powmod48(A_VAL, 1ULL << i);
        INV_POWS_2I[i] = invmod48(p);
    }
}

static u64 dlog_pow2(u64 h) {
    u64 x = 0;
    u64 cur = h & MASK48;
    for (int i = 0; i < ORDER_EXP; i++) {
        u64 e = 1ULL << (ORDER_EXP - 1 - i);
        u64 t = powmod48(cur, e);
        if (t == G_ORDER2) {
            x |= (1ULL << i);
            cur = mul48(cur, INV_POWS_2I[i]);
        }
        // else t should be 1
    }
    return x;
}

// ---- powA_sumY: compute (A^n mod 2^48, Y_n) ----

static void powA_sumY(u64 n, u64 *out_powv, u64 *out_sumv) {
    if (n == 0) {
        *out_powv = 1;
        *out_sumv = 0;
        return;
    }
    u64 powv = 1, sumv = 0;
    // Find highest bit
    int top = 63;
    while (((n >> top) & 1) == 0) top--;
    for (int bit = top; bit >= 0; bit--) {
        // double
        sumv = mul48(sumv, (1 + powv) & MASK48);
        powv = mul48(powv, powv);
        if ((n >> bit) & 1) {
            sumv = (sumv + powv) & MASK48;
            powv = mul48(powv, A_VAL);
        }
    }
    *out_powv = powv;
    *out_sumv = sumv;
}

// ---- index of state ----

static u64 index_of_state(u64 a0, u64 target) {
    a0 &= MASK48;
    target &= MASK48;
    if (target == a0) return 0;

    u64 a1 = step48(a0);
    u64 K = (a1 - a0) & MASK48;
    u64 invK = invmod48(K);

    u64 y_target = mul48((target - a0) & MASK48, invK);
    u64 h = (mul48(A_VAL - 1, y_target) + 1) & MASK48;

    u64 n0 = dlog_pow2(h);
    u64 stepN = 1ULL << ORDER_EXP;

    for (int t = 0; t < 4; t++) {
        u64 n = n0 + (u64)t * stepN;
        u64 pv, sv;
        powA_sumY(n, &pv, &sv);
        if (sv == y_target) return n;
    }
    return 0;  // should not happen
}

// ---- main ----

long long p803_native(void) {
    precompute_dlog();

    // Find unique seed for "PuzzleOne"
    u64 puzzle_states[16];
    int n_puzzle = solve_states("PuzzleOne", puzzle_states, 16);
    // Expect exactly 1
    u64 seed = puzzle_states[0];

    // Find all seeds for "LuckyText"
    u64 lucky_states[256];
    int n_lucky = solve_states("LuckyText", lucky_states, 256);

    u64 best = 0;
    int found = 0;
    for (int i = 0; i < n_lucky; i++) {
        u64 n = index_of_state(seed, lucky_states[i]);
        if (!found || n < best) {
            best = n;
            found = 1;
        }
    }

    return (long long)best;
}
