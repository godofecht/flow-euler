#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

/* Project Euler 967: B-trivisible F(N, B) with N=10^18, B=120.
 *
 * Meet-in-the-middle over the squarefree products of primes <= B (excluding 3).
 * For each squarefree d we track a sign t = (-1)^k and a mod-3 character vector
 * (s0, s1, s2). Combining left and right halves gives the contribution
 *   t_L * t_R * floor(N / (d_L * d_R)) * (s0_L*s0_R + s1_L*s2_R + s2_L*s1_R).
 */

typedef struct {
    u64 prod;
    i64 t;
    i64 s0, s1, s2;
} Item;

static int primes_upto(int B, int *out) {
    char *sieve = calloc(B + 1, 1);
    memset(sieve, 1, B + 1);
    sieve[0] = 0; sieve[1] = 0;
    for (int i = 2; (i64)i * i <= B; i++)
        if (sieve[i])
            for (int j = i * i; j <= B; j += i) sieve[j] = 0;
    int k = 0;
    for (int i = 2; i <= B; i++)
        if (sieve[i] && i != 3) out[k++] = i;
    free(sieve);
    return k;
}

/* Extend the mod-3 character vector when multiplying by prime p (p != 3).
 * cls = p % 3, which is 1 or 2. */
static inline void extend_vec(int cls, i64 s0, i64 s1, i64 s2,
                              i64 *ns0, i64 *ns1, i64 *ns2) {
    if (cls == 1) {
        *ns0 = s0 - s2;
        *ns1 = s1 - s0;
        *ns2 = s2 - s1;
    } else {
        *ns0 = s0 - s1;
        *ns1 = s1 - s2;
        *ns2 = s2 - s0;
    }
}

/* Enumerate squarefree products <= N from primes[idx..) via iterative DFS.
 * Returns array (malloc'd) and count. */
static Item *gen_half(const int *primes, int cnt, u64 N, i64 *out_n) {
    i64 cap = 1 << 20;
    Item *out = malloc(cap * sizeof(Item));
    i64 n = 0;
    /* stack: index, prod, t, s0, s1, s2 */
    i64 scap = 1 << 16;
    i64 *stk = malloc(scap * 6 * sizeof(i64));
    i64 sp = 0;
    stk[sp*6+0] = 0; stk[sp*6+1] = 1; stk[sp*6+2] = 1;
    stk[sp*6+3] = 1; stk[sp*6+4] = 0; stk[sp*6+5] = 0;
    sp++;

    while (sp > 0) {
        sp--;
        i64 idx = stk[sp*6+0];
        u64 prod = (u64)stk[sp*6+1];
        i64 t  = stk[sp*6+2];
        i64 s0 = stk[sp*6+3];
        i64 s1 = stk[sp*6+4];
        i64 s2 = stk[sp*6+5];

        if (n >= cap) {
            cap *= 2;
            out = realloc(out, cap * sizeof(Item));
        }
        out[n].prod = prod;
        out[n].t = t;
        out[n].s0 = s0;
        out[n].s1 = s1;
        out[n].s2 = s2;
        n++;

        for (int j = (int)idx; j < cnt; j++) {
            int p = primes[j];
            if (prod > N / (u64)p) break; /* overflow-safe: prod*p > N */
            u64 np = prod * (u64)p;
            i64 ns0, ns1, ns2;
            extend_vec(p % 3, s0, s1, s2, &ns0, &ns1, &ns2);
            if (sp + 1 > scap) {
                scap *= 2;
                stk = realloc(stk, scap * 6 * sizeof(i64));
            }
            stk[sp*6+0] = j + 1;
            stk[sp*6+1] = (i64)np;
            stk[sp*6+2] = -t;
            stk[sp*6+3] = ns0;
            stk[sp*6+4] = ns1;
            stk[sp*6+5] = ns2;
            sp++;
        }
    }
    free(stk);
    *out_n = n;
    return out;
}

static int cmp_item(const void *a, const void *b) {
    u64 pa = ((const Item *)a)->prod;
    u64 pb = ((const Item *)b)->prod;
    if (pa < pb) return -1;
    if (pa > pb) return 1;
    return 0;
}

/* upper_bound: first index in [lo, hi) with key > val, i.e. bisect_right */
static i64 bisect_right(const u64 *arr, i64 lo, i64 hi, u64 val) {
    while (lo < hi) {
        i64 mid = lo + (hi - lo) / 2;
        if (arr[mid] <= val) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

long long p967_native(void) {
    const u64 N = (u64)1000000000000000000ULL; /* 10^18 */
    const int B = 120;

    int primes[64];
    int cnt = primes_upto(B, primes);

    int mid = cnt / 2;
    const int *L = primes;
    int lcnt = mid;
    const int *R = primes + mid;
    int rcnt = cnt - mid;

    i64 an, bn;
    Item *A = gen_half(L, lcnt, N, &an);
    Item *Bv = gen_half(R, rcnt, N, &bn);

    qsort(Bv, bn, sizeof(Item), cmp_item);

    u64 *bp = malloc(bn * sizeof(u64));
    for (i64 i = 0; i < bn; i++) bp[i] = Bv[i].prod;

    /* prefix sums of t*s0, t*s1, t*s2 */
    i64 *pref0 = malloc((bn + 1) * sizeof(i64));
    i64 *pref1 = malloc((bn + 1) * sizeof(i64));
    i64 *pref2 = malloc((bn + 1) * sizeof(i64));
    pref0[0] = pref1[0] = pref2[0] = 0;
    for (i64 i = 0; i < bn; i++) {
        i64 t = Bv[i].t;
        pref0[i+1] = pref0[i] + t * Bv[i].s0;
        pref1[i+1] = pref1[i] + t * Bv[i].s1;
        pref2[i+1] = pref2[i] + t * Bv[i].s2;
    }

    i128 ans = 0;
    for (i64 a = 0; a < an; a++) {
        u64 ap = A[a].prod;
        if (ap > N) continue;
        u64 limit = N / ap;
        i64 r = bisect_right(bp, 0, bn, limit);
        i64 i = 0;
        while (i < r) {
            u64 v = limit / bp[i];
            if (v == 0) break;
            u64 maxb = limit / v;
            i64 j = bisect_right(bp, i, r, maxb);
            i64 s0 = pref0[j] - pref0[i];
            i64 s1 = pref1[j] - pref1[i];
            i64 s2 = pref2[j] - pref2[i];
            i128 comb = (i128)A[a].s0 * s0 + (i128)A[a].s1 * s2 + (i128)A[a].s2 * s1;
            ans += (i128)A[a].t * (i64)v * comb;
            i = j;
        }
    }

    free(A); free(Bv); free(bp);
    free(pref0); free(pref1); free(pref2);
    return (long long)ans;
}
