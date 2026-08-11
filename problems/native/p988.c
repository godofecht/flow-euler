#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef int64_t i64;

/* Project Euler 988
 * frog_sum(19, 53) = 2727531976556215755
 *
 * DP keyed by height (0..~48). Each entry holds (count, total).
 * Implemented with flat arrays since the height range is small.
 */

#define MAXW 64
#define MAXH 64

i64 p988_native(void) {
    i64 a = 19, b = 53;
    if (a > b) { i64 t = a; a = b; b = t; }
    if (a == 1) return 0;

    i64 width = b - 1;
    i64 h[MAXW + 1];
    for (i64 i = 1; i <= width; i++) {
        h[i] = (a * b - a * i - 1) / b;
    }

    /* dp arrays indexed by height 0..MAXH */
    i64 cnt[MAXH + 1];
    i64 tot[MAXH + 1];
    i64 ncnt[MAXH + 1];
    i64 ntot[MAXH + 1];

    memset(cnt, 0, sizeof(cnt));
    memset(tot, 0, sizeof(tot));

    for (i64 t = 0; t <= h[1]; t++) {
        cnt[t] = 1;
        tot[t] = 0;
    }

    for (i64 i = 2; i <= width; i++) {
        memset(ncnt, 0, sizeof(ncnt));
        memset(ntot, 0, sizeof(ntot));
        for (i64 prev = 0; prev <= MAXH; prev++) {
            if (cnt[prev] == 0 && tot[prev] == 0) continue;
            i64 count = cnt[prev];
            i64 total = tot[prev];
            i64 limit = prev < h[i] ? prev : h[i];
            for (i64 cur = 0; cur <= limit; cur++) {
                i64 add = 0;
                if (prev > cur && prev > 0) {
                    add = a * b - a * (i - 1) - b * prev;
                }
                ncnt[cur] += count;
                ntot[cur] += total + count * add;
            }
        }
        memcpy(cnt, ncnt, sizeof(cnt));
        memcpy(tot, ntot, sizeof(tot));
    }

    i64 answer = 0;
    i64 last_column = width;
    for (i64 prev = 0; prev <= MAXH; prev++) {
        if (cnt[prev] == 0 && tot[prev] == 0) continue;
        i64 add = 0;
        if (prev > 0) {
            add = a * b - a * last_column - b * prev;
        }
        answer += tot[prev] + cnt[prev] * add;
    }
    return answer;
}
