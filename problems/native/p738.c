#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef int64_t i64;
typedef __int128 i128;

static const i64 MOD = 1000000007LL;

/* Problem 738: Counting Ordered Factorisations.
   D(10^10, 10^10) mod 1e9+7.
   Port of the Python reference solver. */

static i64 icbrt_floor(i64 n) {
    if (n <= 0) return 0;
    double x = cbrt((double)n);
    i64 r = (i64)x;
    /* correct */
    while ((r + 1) * (r + 1) * (r + 1) <= n) r++;
    while (r * r * r > n) r--;
    return r;
}

static i64 sum_floor_range(i64 m, i64 l, i64 r) {
    i64 res = 0;
    i64 i = l;
    while (i <= r) {
        i64 q = m / i;
        i64 j = m / q;
        if (j > r) j = r;
        res += q * (j - i + 1);
        i = j + 1;
    }
    return res;
}

static i64 sum_arith(i64 l, i64 r) {
    i64 n = r - l + 1;
    return (l + r) * n / 2;
}

/* Memoization for _count_and_length(m, a).
   The recursion is on (m, a) where a >= 2 and m <= N=1e10.
   We use a hash table keyed by (m, a). */

typedef struct {
    i64 m;
    i64 a;
    i64 c;   /* count mod MOD */
    i64 l;   /* length sum mod MOD */
} Entry;

#define HASH_CAP (1 << 22)
static Entry *htab = NULL;

static i64 hash_key(i64 m, i64 a) {
    /* m up to 1e10, a up to ~1e5; mix */
    i64 h = m * 1000003LL + a;
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdLL;
    h ^= h >> 33;
    return h;
}

static void htab_init(void) {
    if (!htab) {
        htab = calloc(HASH_CAP, sizeof(Entry));
    } else {
        memset(htab, 0, HASH_CAP * sizeof(Entry));
    }
}

/* returns 1 if found, fills *cp,*lp */
static int htab_get(i64 m, i64 a, i64 *cp, i64 *lp) {
    i64 h = hash_key(m, a);
    i64 idx = h & (HASH_CAP - 1);
    while (1) {
        if (htab[idx].m == 0 && htab[idx].a == 0) {
            /* empty (but m=0,a=0 is never a real query since a>=2) */
            return 0;
        }
        if (htab[idx].m == m && htab[idx].a == a) {
            *cp = htab[idx].c;
            *lp = htab[idx].l;
            return 1;
        }
        idx = (idx + 1) & (HASH_CAP - 1);
    }
}

static void htab_put(i64 m, i64 a, i64 c, i64 l) {
    i64 h = hash_key(m, a);
    i64 idx = h & (HASH_CAP - 1);
    while (1) {
        if (htab[idx].m == 0 && htab[idx].a == 0) {
            htab[idx].m = m;
            htab[idx].a = a;
            htab[idx].c = c;
            htab[idx].l = l;
            return;
        }
        if (htab[idx].m == m && htab[idx].a == a) {
            htab[idx].c = c;
            htab[idx].l = l;
            return;
        }
        idx = (idx + 1) & (HASH_CAP - 1);
    }
}

static void count_and_length(i64 m, i64 a, i64 *Cp, i64 *Lp);

static void count_and_length(i64 m, i64 a, i64 *Cp, i64 *Lp) {
    if (m < a) { *Cp = 0; *Lp = 0; return; }

    i64 Cc, Ll;
    if (htab_get(m, a, &Cc, &Ll)) {
        *Cp = Cc; *Lp = Ll;
        return;
    }

    i64 aa = a * a;
    i64 C, L;
    if (aa > m) {
        i64 cnt = (m - a + 1) % MOD;
        C = cnt; L = cnt;
        goto done;
    }

    /* length-1 tuples */
    C = m - a + 1;
    L = C;

    i64 s = (i64)sqrtl((long double)m);
    /* fix sqrt */
    while ((s + 1) * (s + 1) <= m) s++;
    while (s * s > m) s--;

    if (aa * a > m) {
        i64 l = a;
        if (l <= s) {
            i64 sf = sum_floor_range(m, l, s);
            i64 sa = sum_arith(l, s);
            i64 baseC = sf - sa + (s - l + 1);
            C += baseC;
            L += 2 * baseC;
        }
        goto done;
    }

    i64 t = icbrt_floor(m);
    i64 upto = (t < s) ? t : s;

    /* recursive part: f in [a..upto] */
    if (upto >= a) {
        for (i64 f = a; f <= upto; f++) {
            i64 subC, subL;
            count_and_length(m / f, f, &subC, &subL);
            C += subC;
            L += subL + subC;
        }
    }

    /* base part: f in [max(a, upto+1)..s] */
    {
        i64 l = (a > upto + 1) ? a : (upto + 1);
        if (l <= s) {
            i64 sf = sum_floor_range(m, l, s);
            i64 sa = sum_arith(l, s);
            i64 baseC = sf - sa + (s - l + 1);
            C += baseC;
            L += 2 * baseC;
        }
    }

done:
    C %= MOD; L %= MOD;
    htab_put(m, a, C, L);
    *Cp = C; *Lp = L;
}

long long p738_native(void) {
    htab_init();
    i64 N = 10000000000LL; /* 10^10 */
    i64 K = 10000000000LL;
    i64 C, L;
    count_and_length(N, 2, &C, &L);
    i64 ans = ((K % MOD) + ((K + 1) % MOD) * C % MOD - L) % MOD;
    ans %= MOD;
    if (ans < 0) ans += MOD;
    return (long long)ans;
}
