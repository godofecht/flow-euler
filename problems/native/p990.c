// Project Euler 990: Digit Equations
// DP over decimal columns with signed carry, mod 1e9+7.
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef long long i64;
typedef __int128 i128;

#define MOD 1000000007LL
#define MAX_N 50
#define MAX_TERMS ((MAX_N + 1) / 2)  // 25
#define MAX_CARRY 25
#define CARRY_RANGE 51  // -25..25

static i64 mod_add(i64 a, i64 b) {
    a += b;
    if (a >= MOD) a -= MOD;
    return a;
}

static i64 mod_mul(i64 a, i64 b) {
    return (i64)((i128)a * b % MOD);
}

// Binomial coefficients C[n][k] mod MOD, for 0 <= n,k <= MAX_TERMS.
static i64 binom[26][26];

static void build_binom() {
    for (int n = 0; n <= MAX_TERMS; n++) {
        binom[n][0] = binom[n][n] = 1;
        for (int k = 1; k < n; k++)
            binom[n][k] = mod_add(binom[n-1][k-1], binom[n-1][k]);
    }
}

// Sum tables: F[p][q][m] = ways to write m as sum of p vars in [0,9] and q vars in [0,8].
// Max degree = 9*p + 8*q. For p,q <= 2*MAX_TERMS=50, max degree = 850.
#define MAX_DEG 851
static i64 sum_tables[51][51][MAX_DEG];
static int sum_table_len[51][51];

static void convolve_small(const i64 *poly, int plen, int width, i64 *out, int *out_len) {
    int n = plen + width - 1;
    for (int i = 0; i < n; i++) out[i] = 0;
    for (int i = 0; i < plen; i++) {
        if (poly[i] == 0) continue;
        for (int d = 0; d < width; d++)
            out[i + d] = mod_add(out[i + d], poly[i]);
    }
    *out_len = n;
}

static void build_sum_tables() {
    // ways_0_to_9[p] = convolution of p copies of [1]*10
    static i64 ways[51][MAX_DEG];
    static i64 tmp[MAX_DEG];
    int ways_len[51];

    ways[0][0] = 1; ways_len[0] = 1;
    for (int p = 1; p <= 2 * MAX_TERMS; p++) {
        convolve_small(ways[p-1], ways_len[p-1], 10, ways[p], &ways_len[p]);
    }

    for (int p = 0; p <= 2 * MAX_TERMS; p++) {
        memcpy(sum_tables[p][0], ways[p], ways_len[p] * sizeof(i64));
        sum_table_len[p][0] = ways_len[p];
        for (int q = 1; q <= 2 * MAX_TERMS; q++) {
            convolve_small(sum_tables[p][q-1], sum_table_len[p][q-1], 9,
                          sum_tables[p][q], &sum_table_len[p][q]);
        }
    }
}

// Transition: (next_left, next_right, next_carry, weight)
typedef struct {
    int nl, nr, nc;
    i64 w;
} Trans;

static Trans *trans_mem[26][26][CARRY_RANGE];
static int trans_count[26][26][CARRY_RANGE];
static int trans_computed[26][26][CARRY_RANGE];

static void compute_transitions(int al, int ar, int carry) {
    int ci = carry + MAX_CARRY;
    if (trans_computed[al][ar][ci]) return;
    trans_computed[al][ar][ci] = 1;

    if (al == 0 && ar == 0) {
        trans_count[al][ar][ci] = 0;
        return;
    }

    // First pass: count transitions.
    int count = 0;
    for (int nl = 0; nl <= al; nl++) {
        for (int nr = 0; nr <= ar; nr++) {
            int continuing = nl + nr;
            int ending = (al - nl) + (ar - nr);
            int ending_left = al - nl;
            const i64 *counts = sum_tables[continuing][ending];
            int clen = sum_table_len[continuing][ending];
            i64 base = -carry - ending_left + 9 * ar;
            for (int nc = -MAX_CARRY; nc <= MAX_CARRY; nc++) {
                int index = 10 * nc + (int)base;
                if (index >= 0 && index < clen && counts[index]) {
                    count++;
                }
            }
        }
    }

    Trans *list = malloc(count * sizeof(Trans));
    int idx = 0;
    for (int nl = 0; nl <= al; nl++) {
        i64 choose_left = binom[al][nl];
        int ending_left = al - nl;
        for (int nr = 0; nr <= ar; nr++) {
            i64 choose_terms = mod_mul(choose_left, binom[ar][nr]);
            int continuing = nl + nr;
            int ending = ending_left + (ar - nr);
            const i64 *counts = sum_tables[continuing][ending];
            int clen = sum_table_len[continuing][ending];
            i64 base = -carry - ending_left + 9 * ar;
            for (int nc = -MAX_CARRY; nc <= MAX_CARRY; nc++) {
                int index = 10 * nc + (int)base;
                if (index >= 0 && index < clen && counts[index]) {
                    i64 weight = mod_mul(choose_terms, counts[index]);
                    list[idx].nl = nl;
                    list[idx].nr = nr;
                    list[idx].nc = nc;
                    list[idx].w = weight;
                    idx++;
                }
            }
        }
    }

    trans_mem[al][ar][ci] = list;
    trans_count[al][ar][ci] = count;
}

// DP: dp[length][al][ar][carry+25]
static i64 dp[51][26][26][CARRY_RANGE];

static i64 solve(int limit) {
    memset(dp, 0, sizeof(dp));

    // Initial states: choose number of terms on each side.
    for (int lt = 1; lt <= MAX_TERMS; lt++) {
        for (int rt = 1; rt <= MAX_TERMS; rt++) {
            int base_length = lt + rt - 1;
            if (base_length <= limit) {
                dp[base_length][lt][rt][MAX_CARRY] = mod_add(
                    dp[base_length][lt][rt][MAX_CARRY], 1);
            }
        }
    }

    i64 answer = 0;
    for (int used = 0; used <= limit; used++) {
        // Add (0,0,0) states to answer.
        answer = mod_add(answer, dp[used][0][0][MAX_CARRY]);

        for (int al = 0; al <= MAX_TERMS; al++) {
            for (int ar = 0; ar <= MAX_TERMS; ar++) {
                if (al == 0 && ar == 0) continue;
                for (int ci = 0; ci < CARRY_RANGE; ci++) {
                    i64 ways = dp[used][al][ar][ci];
                    if (ways == 0) continue;
                    int carry = ci - MAX_CARRY;

                    compute_transitions(al, ar, carry);
                    Trans *list = trans_mem[al][ar][ci];
                    int cnt = trans_count[al][ar][ci];

                    int next_length = used + al + ar;
                    if (next_length > limit) continue;

                    for (int t = 0; t < cnt; t++) {
                        int nci = list[t].nc + MAX_CARRY;
                        i64 val = mod_mul(ways, list[t].w);
                        dp[next_length][list[t].nl][list[t].nr][nci] =
                            mod_add(dp[next_length][list[t].nl][list[t].nr][nci], val);
                    }
                }
            }
        }
    }

    return answer;
}

long long p990_native(void) {
    build_binom();
    build_sum_tables();
    return solve(MAX_N);
}
