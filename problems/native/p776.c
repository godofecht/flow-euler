/* Project Euler 776: Digit Sum Division
 *
 * For n >= 1, let d(n) be the sum of decimal digits of n.
 * F(N) = sum_{n=1..N} n / d(n).
 *
 * A digit DP over the decimal representation of N groups every integer in
 * [0, N] by its digit sum s.  For each s we track the count of integers and
 * their sum.  F(N) = sum_{s>=1} (sum of integers with digit sum s) / s.
 *
 * The integer sums reach ~1e36, beyond exact double range, so the DP carries
 * the sums in long double (64-bit mantissa on x86) and only the final
 * quotient is narrowed to double.
 */

#include <stdio.h>
#include <string.h>

double p776_native(void) {
    /* N = 1234567890123456789 */
    static const int digits[] = {
        1, 2, 3, 4, 5, 6, 7, 8, 9,
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9
    };
    const int L = (int)(sizeof(digits) / sizeof(digits[0]));
    const int max_sum = 9 * L; /* 171 */

    /* tight: prefixes still equal to N; loose: already strictly smaller. */
    long double cnt_tight[172], sum_tight[172];
    long double cnt_loose[172], sum_loose[172];
    long double ncnt_tight[172], nsum_tight[172];
    long double ncnt_loose[172], nsum_loose[172];

    memset(cnt_tight, 0, sizeof(cnt_tight));
    memset(sum_tight, 0, sizeof(sum_tight));
    memset(cnt_loose, 0, sizeof(cnt_loose));
    memset(sum_loose, 0, sizeof(sum_loose));
    cnt_tight[0] = 1.0L;

    for (int i = 0; i < L; i++) {
        int lim = digits[i];

        memset(ncnt_tight, 0, sizeof(ncnt_tight));
        memset(nsum_tight, 0, sizeof(nsum_tight));
        memset(ncnt_loose, 0, sizeof(ncnt_loose));
        memset(nsum_loose, 0, sizeof(nsum_loose));

        /* Loose -> loose: next digit free (0..9). */
        for (int s = 0; s <= max_sum; s++) {
            long double c = cnt_loose[s];
            if (c == 0.0L) continue;
            long double v10 = sum_loose[s] * 10.0L;
            for (int d = 0; d <= 9; d++) {
                int ns = s + d;
                if (ns > max_sum) break;
                ncnt_loose[ns] += c;
                nsum_loose[ns] += v10 + c * (long double)d;
            }
        }

        /* Tight -> tight/loose: next digit restricted by lim. */
        for (int s = 0; s <= max_sum; s++) {
            long double c = cnt_tight[s];
            if (c == 0.0L) continue;
            long double v10 = sum_tight[s] * 10.0L;
            for (int d = 0; d <= lim; d++) {
                int ns = s + d;
                if (ns > max_sum) break;
                if (d == lim) {
                    ncnt_tight[ns] += c;
                    nsum_tight[ns] += v10 + c * (long double)d;
                } else {
                    ncnt_loose[ns] += c;
                    nsum_loose[ns] += v10 + c * (long double)d;
                }
            }
        }

        memcpy(cnt_tight, ncnt_tight, sizeof(cnt_tight));
        memcpy(sum_tight, nsum_tight, sizeof(sum_tight));
        memcpy(cnt_loose, ncnt_loose, sizeof(cnt_loose));
        memcpy(sum_loose, nsum_loose, sizeof(sum_loose));
    }

    /* F(N) = sum_{s>=1} (sum_tight[s] + sum_loose[s]) / s.
     * s = 0 only matches the number 0, which contributes nothing. */
    long double total = 0.0L;
    for (int s = 1; s <= max_sum; s++) {
        long double sums = sum_tight[s] + sum_loose[s];
        if (sums != 0.0L) {
            total += sums / (long double)s;
        }
    }

    return (double)total;
}

void p776_print(void) {
    double v = p776_native();
    char buf[64];
    snprintf(buf, sizeof(buf), "%.12e\n", v);
    /* Remove the '+' in exponent to match expected format e.g. e33 not e+33 */
    char *e = strchr(buf, 'e');
    if (e && e[1] == '+') {
        memmove(e + 1, e + 2, strlen(e + 2) + 1);
    }
    fputs(buf, stdout);
}
