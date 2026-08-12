// Project Euler 867: Tiling Dodecagon
// T(10) mod 1e9+7
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
typedef __int128 i128;

#define MOD 1000000007LL

#define MAX_L 19  // 2*10-1

static i64 indep_count[MAX_L + 1];
static i64 indep_lists[MAX_L + 1][1 << MAX_L]; // over-allocated; we index by count

static i64 modpow(i64 base, i64 exp, i64 mod) {
    i64 r = 1; base %= mod;
    while (exp > 0) {
        if (exp & 1) r = (i64)((i128)r * base % mod);
        base = (i64)((i128)base * base % mod);
        exp >>= 1;
    }
    return r;
}

static void init_indep(void) {
    for (int L = 0; L <= MAX_L; L++) {
        i64 cnt = 0;
        for (i64 m = 0; m < (1LL << L); m++) {
            if ((m & (m << 1)) == 0)
                indep_lists[L][cnt++] = m;
        }
        indep_count[L] = cnt;
    }
}

// Subset-sum (zeta) transform: out[mask] = sum(values[sub] for sub subset mask) mod MOD
// Runs in O(L * 2^L)
static void subset_sums(i64 *out, const i64 *values, int L) {
    i64 size = 1LL << L;
    memcpy(out, values, (size_t)size * sizeof(i64));
    for (int i = 0; i < L; i++) {
        i64 step = 1LL << i;
        i64 block = step << 1;
        for (i64 start = 0; start < size; start += block) {
            i64 mid = start + step;
            i64 end = start + block;
            for (i64 m = mid; m < end; m++) {
                i64 s = out[m] + out[m - step];
                out[m] = (s >= MOD) ? s - MOD : s;
            }
        }
    }
}

// Count independent sets on a layered lattice described by row lengths.
// row_lengths: array of ints, each consecutive pair differs by exactly 1.
static i64 count_independent_sets(const int *row_lengths, int nrows) {
    if (nrows == 0) return 1;

    int L0 = row_lengths[0];
    i64 size0 = 1LL << L0;
    i64 *dp = (i64*)calloc((size_t)size0, sizeof(i64));
    for (i64 c = 0; c < indep_count[L0]; c++)
        dp[indep_lists[L0][c]] = 1;

    for (int i = 1; i < nrows; i++) {
        int Lc = row_lengths[i - 1];
        int Ln = row_lengths[i];
        i64 size_c = 1LL << Lc;
        i64 size_n = 1LL << Ln;

        i64 *subs = (i64*)malloc((size_t)size_c * sizeof(i64));
        subset_sums(subs, dp, Lc);

        i64 *dp2 = (i64*)calloc((size_t)size_n, sizeof(i64));
        i64 fullmask = size_c - 1;

        if (Ln == Lc + 1) {
            for (i64 c = 0; c < indep_count[Ln]; c++) {
                i64 b = indep_lists[Ln][c];
                i64 forb = (b | (b >> 1)) & fullmask;
                i64 allowed = fullmask ^ forb;
                dp2[b] = subs[allowed];
            }
        } else { // Ln == Lc - 1
            for (i64 c = 0; c < indep_count[Ln]; c++) {
                i64 b = indep_lists[Ln][c];
                i64 forb = (b | (b << 1)) & fullmask;
                i64 allowed = fullmask ^ forb;
                dp2[b] = subs[allowed];
            }
        }

        free(subs);
        free(dp);
        dp = dp2;
    }

    int last_L = row_lengths[nrows - 1];
    i64 total = 0;
    for (i64 c = 0; c < indep_count[last_L]; c++) {
        total += dp[indep_lists[last_L][c]];
        if (total >= MOD) total -= MOD;
    }
    free(dp);
    return total % MOD;
}

// H(n): tilings of regular hexagon of side n
static i64 H(int n) {
    int inc[32], dec[32];
    int ni = 0, nd = 0;
    for (int i = n; i < 2 * n; i++) inc[ni++] = i;
    for (int i = 2 * n - 2; i >= n; i--) dec[nd++] = i;
    int rows[64];
    memcpy(rows, inc, ni * sizeof(int));
    memcpy(rows + ni, dec, nd * sizeof(int));
    return count_independent_sets(rows, ni + nd);
}

// F(n, h): tilings of truncated equilateral triangle
static i64 F(int n, int h) {
    int rows = h - 1;
    if (rows <= 0) return 1;
    int lengths[64];
    for (int i = 0; i < rows; i++) {
        int v = (n - 2) - i;
        lengths[i] = v > 0 ? v : 0;
    }
    // Remove trailing zeros (count_independent_sets handles 0-length rows as return 1)
    // Actually we need to handle this: if any length is 0, the lattice has no nodes there
    // But the Python code passes max(0, ...) and count_independent_sets returns 1 for empty
    // Let's filter: find first 0 and truncate
    int eff_rows = rows;
    for (int i = 0; i < rows; i++) {
        if (lengths[i] == 0) {
            eff_rows = i;
            break;
        }
    }
    return count_independent_sets(lengths, eff_rows);
}

// R(u, v): helper recursion with memoization
// R values indexed by (u,v) with u,v <= 10
static i64 R_cache[16][16];
static int R_computed[16][16];

static i64 R(int u, int v);

static i64 R(int u, int v) {
    if (v == 0) return H(u);
    if (R_computed[u][v]) return R_cache[u][v];

    i64 res = (u == 1 && v == 1) ? 1 : 0;
    for (int w = 0; w < u; w++) {
        i64 corner = F(u, u - w);
        i64 r_vw = R(v, w);
        i64 c6 = modpow(corner, 6, MOD);
        res = (res + (i64)((i128)r_vw * c6 % MOD)) % MOD;
    }
    R_computed[u][v] = 1;
    R_cache[u][v] = res;
    return res;
}

// T(n): tilings of regular dodecagon of side n
static i64 T(int n) {
    i64 ans = (2 * R(n, n) - (n == 1 ? 1 : 0)) % MOD;
    if (ans < 0) ans += MOD;
    return ans;
}

long long p867_native(void) {
    init_indep();
    return T(10);
}
