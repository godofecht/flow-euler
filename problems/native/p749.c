#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t  i64;
typedef __int128 i128;

/* Problem 749: Near Power Sums.
   Sum of all near power sum numbers with at most 16 digits.
   Port of the Python reference solver. */

#define MAX_DIGITS 16
#define POWER_CAP ((i128)1000000000000000000LL) /* 10^18, cap to avoid overflow */

static i64 p10[MAX_DIGITS + 1];
static int max_k;
static i128 pow_table[64][10];   /* pow_table[k][d] = d^k, capped at POWER_CAP */

/* k bounds */
static int k_low[MAX_DIGITS + 1][10];
static int k_high[MAX_DIGITS + 1][10][MAX_DIGITS + 1];

/* digit packing */
static i64 pack4[10000];
static i64 pack_exact[5][10000]; /* pack_exact[len][x], len 1..4 */

/* ---- hash set for results ---- */
#define SET_CAP (1 << 16)
static i64 *set_tbl = NULL;

static void set_init(void) {
    if (!set_tbl) set_tbl = malloc(SET_CAP * sizeof(i64));
    memset(set_tbl, 0, SET_CAP * sizeof(i64));
}

static int set_insert(i64 v) {
    /* 0 is never a result (results are >= 10^(L-1) >= 1) */
    i64 h = v;
    h ^= h >> 33; h *= 0xff51afd7ed558ccdLL; h ^= h >> 33;
    i64 idx = h & (SET_CAP - 1);
    while (1) {
        if (set_tbl[idx] == 0) { set_tbl[idx] = v; return 1; }
        if (set_tbl[idx] == v) return 0;
        idx = (idx + 1) & (SET_CAP - 1);
    }
}

/* ---- build tables ---- */
static void build_pow10(void) {
    p10[0] = 1;
    for (int i = 1; i <= MAX_DIGITS; i++) p10[i] = p10[i-1] * 10;
}

static void compute_max_k(void) {
    i64 limit = p10[MAX_DIGITS];
    int k = 0; i64 p = 1;
    while (p <= limit) { k++; p *= 2; }
    max_k = k;
}

static void build_pow_table(void) {
    for (int d = 0; d < 10; d++) pow_table[1][d] = d;
    for (int k = 2; k <= max_k; k++) {
        for (int d = 0; d < 10; d++) {
            i128 val = pow_table[k-1][d] * d;
            if (val > POWER_CAP) val = POWER_CAP;
            pow_table[k][d] = val;
        }
    }
}

static void build_k_bounds(void) {
    for (int L = 1; L <= MAX_DIGITS; L++) {
        i64 lt = p10[L-1] - 1;
        i64 ht = p10[L];
        for (int m = 0; m < 10; m++) {
            k_low[L][m] = max_k + 1;
            for (int c = 0; c <= MAX_DIGITS; c++)
                k_high[L][m][c] = 0;
        }
        for (int m = 2; m < 10; m++) {
            /* k_low: smallest k with L*m^k >= lt */
            i128 p = m; int k = 1;
            while (k <= max_k && (i128)L * p < lt) { p *= m; k++; }
            k_low[L][m] = k;

            /* k_high: largest k with c*m^k <= ht */
            for (int c = 1; c <= L; c++) {
                i128 p2 = m; int best = 0;
                for (int kk = 1; kk <= max_k; kk++) {
                    if ((i128)c * p2 <= ht) { best = kk; p2 *= m; }
                    else break;
                }
                k_high[L][m][c] = best;
            }
        }
    }
}

static void build_pack_tables(void) {
    for (int x = 0; x < 10000; x++) {
        i64 code = 0, y = x;
        for (int i = 0; i < 4; i++) {
            int d = y % 10; code += (i64)1 << (5*d); y /= 10;
        }
        pack4[x] = code;
    }
    for (int x = 0; x < 10; x++) pack_exact[1][x] = (i64)1 << (5*x);
    for (int x = 0; x < 100; x++) {
        i64 code = 0, y = x;
        for (int i = 0; i < 2; i++) { int d = y%10; code += (i64)1<<(5*d); y/=10; }
        pack_exact[2][x] = code;
    }
    for (int x = 0; x < 1000; x++) {
        i64 code = 0, y = x;
        for (int i = 0; i < 3; i++) { int d = y%10; code += (i64)1<<(5*d); y/=10; }
        pack_exact[3][x] = code;
    }
    for (int x = 0; x < 10000; x++) pack_exact[4][x] = pack4[x];
}

static i64 pack_digits_len(i64 n, int L) {
    if (L <= 4) return pack_exact[L][n];
    if (L <= 8) {
        i64 a = n % 10000; i64 b = n / 10000;
        return pack4[a] + pack_exact[L-4][b];
    }
    if (L <= 12) {
        i64 a = n % 10000; n /= 10000;
        i64 b = n % 10000; i64 c = n / 10000;
        return pack4[a] + pack4[b] + pack_exact[L-8][c];
    }
    /* L <= 16 */
    i64 a = n % 10000; n /= 10000;
    i64 b = n % 10000; n /= 10000;
    i64 c = n % 10000; i64 d = n / 10000;
    return pack4[a] + pack4[b] + pack4[c] + pack_exact[L-12][d];
}

/* ---- enumeration state (globals for leaf/rec) ---- */
static int   g_Ks[64];   /* k indices */
static int   g_nks;      /* number of k values */
static i128  g_base_t[64];
static int   g_counts[10];
static i64   g_sig_m;
static i64   g_L_lo, g_L_hi;
static int   g_L;

static void leaf(i64 packed_part) {
    i64 sig = packed_part | g_sig_m;
    for (int idx = 0; idx < g_nks; idx++) {
        i128 t = g_base_t[idx];
        int k = g_Ks[idx];
        for (int d = 0; d < 10; d++) {
            int c = g_counts[d];
            if (c) t += (i128)c * pow_table[k][d];
        }
        /* n = t - 1 */
        if (t > 0) {
            i64 n = (i64)(t - 1);
            if (g_L_lo <= n && n < g_L_hi && pack_digits_len(n, g_L) == sig)
                set_insert(n);
        }
        /* n = t + 1 */
        {
            i64 n = (i64)(t + 1);
            if (g_L_lo <= n && n < g_L_hi && pack_digits_len(n, g_L) == sig)
                set_insert(n);
        }
    }
}

static void rec(int d, int remaining, i64 packed_part, int m) {
    if (d == m - 1) {
        g_counts[d] = remaining;
        leaf(packed_part | ((i64)remaining << (5 * d)));
        return;
    }
    int shift = 5 * d;
    for (int c = 0; c <= remaining; c++) {
        g_counts[d] = c;
        rec(d + 1, remaining - c, packed_part | ((i64)c << shift), m);
    }
}

long long p749_native(void) {
    build_pow10();
    compute_max_k();
    build_pow_table();
    build_k_bounds();
    build_pack_tables();
    set_init();

    for (int L = 1; L <= MAX_DIGITS; L++) {
        g_L = L;
        g_L_lo = p10[L-1];
        g_L_hi = p10[L];

        for (int m = 2; m < 10; m++) {
            int shift_m = 5 * m;
            int km_base = k_low[L][m];

            for (int c_m = 1; c_m <= L; c_m++) {
                int k1 = km_base;
                int k2 = k_high[L][m][c_m];
                if (k1 > k2) continue;

                g_nks = 0;
                for (int k = k1; k <= k2; k++) {
                    g_Ks[g_nks] = k;
                    g_base_t[g_nks] = (i128)c_m * pow_table[k][m];
                    g_nks++;
                }
                int rem = L - c_m;
                for (int i = 0; i < m; i++) g_counts[i] = 0;
                g_sig_m = (i64)c_m << shift_m;

                if (m == 2) {
                    for (int c0 = 0; c0 <= rem; c0++) {
                        int c1 = rem - c0;
                        g_counts[0] = c0; g_counts[1] = c1;
                        leaf((i64)c0 | ((i64)c1 << 5));
                    }
                } else {
                    rec(0, rem, 0, m);
                }
            }
        }
    }

    /* sum all results */
    i128 total = 0;
    for (int i = 0; i < SET_CAP; i++) {
        if (set_tbl[i]) total += set_tbl[i];
    }
    return (long long)total;
}
