// Project Euler 974
// Count numbers using digits {1,3,5,7,9}, divisible by 21, using all five
// digits, ending in 5. Find the 10^16-th such number (ordered by length then
// lexicographically). The answer is a 29-digit string that exceeds i64, so
// the native helper prints it directly and returns 0.
#include <stdint.h>
#include <stdio.h>
#include <string.h>

typedef __int128 i128;

static const int DIGITS[5] = {1, 3, 5, 7, 9};
/* BIT[d] maps digit d to bit index 0..4 */
static const int BIT[10] = {-1, 0, -1, 1, -1, 2, -1, 3, -1, 4};
#define ALL 31
#define MAXL 201

/* dp[pos][mod3][mod7][mask]: number of valid completions from pos to L */
static i128 dp[MAXL][3][7][32];

static i128 count_len(int L) {
    memset(dp, 0, sizeof(dp));

    /* Base case: at position L, valid only if mod3==0, mod7==0, mask==ALL */
    dp[L][0][0][ALL] = (i128)1;

    for (int pos = L - 1; pos >= 0; pos--) {
        int nchoices = (pos == L - 1) ? 1 : 5;
        for (int m3 = 0; m3 < 3; m3++) {
            for (int m7 = 0; m7 < 7; m7++) {
                for (int mask = 0; mask < 32; mask++) {
                    i128 total = 0;
                    for (int ci = 0; ci < nchoices; ci++) {
                        int d = (pos == L - 1) ? 5 : DIGITS[ci];
                        int nm3 = (m3 * 10 + d) % 3;
                        int nm7 = (m7 * 10 + d) % 7;
                        int nmask = mask ^ (1 << BIT[d]);
                        total += dp[pos + 1][nm3][nm7][nmask];
                    }
                    dp[pos][m3][m7][mask] = total;
                }
            }
        }
    }
    return dp[0][0][0][0];
}

static void unrank(int L, i128 k, char *out) {
    int m3 = 0, m7 = 0, mask = 0;
    for (int pos = 0; pos < L; pos++) {
        int nchoices = (pos == L - 1) ? 1 : 5;
        for (int ci = 0; ci < nchoices; ci++) {
            int d = (pos == L - 1) ? 5 : DIGITS[ci];
            int nm3 = (m3 * 10 + d) % 3;
            int nm7 = (m7 * 10 + d) % 7;
            int nmask = mask ^ (1 << BIT[d]);
            i128 cnt = dp[pos + 1][nm3][nm7][nmask];
            if (k > cnt) {
                k -= cnt;
            } else {
                out[pos] = (char)('0' + d);
                m3 = nm3;
                m7 = nm7;
                mask = nmask;
                break;
            }
        }
    }
    out[L] = '\0';
}

long long p974_native(void) {
    /* n = 10^16 */
    i128 n = (i128)1;
    for (int i = 0; i < 16; i++) n *= 10;

    i128 cum = 0;
    int targetL = -1;
    for (int L = 1; L <= 200; L += 2) {
        i128 c = count_len(L);
        if (cum + c >= n) {
            targetL = L;
            break;
        }
        cum += c;
    }

    /* dp is now filled for targetL */
    char out[256];
    unrank(targetL, n - cum, out);
    printf("%s\n", out);
    return 0;
}
