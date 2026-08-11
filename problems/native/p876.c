/* Project Euler 876: Triplet Tricks.
   Computes sum_{k=1}^{18} F(6^k, 10^k) where F(a,b) = sum_{c>=1} f(a,b,c).

   f(a,b,c) is the minimum number of involutions until an entry hits 0.
   Candidate c values come from divisor pairs (x|a, y|b):
       c1 = (x+y)(a/x + b/y)
       c2 = (x-y)(a/x - b/y)   (only if positive)
   The step count is the subtractive-Euclidean step count for (x,y),
   minus one for the second family.

   For k=18, a=6^18 and b=10^18, so c can be as large as ~a*b ~ 1e32,
   which overflows int64. The c keys are stored as __int128; the step
   counts and the final sum fit in int64.
*/

#include <stdint.h>
#include <stdlib.h>

typedef long long i64;
typedef __int128 i128;

static i64 euclid_subtraction_steps(i64 m, i64 n) {
    i64 steps = 0;
    while (n) {
        i64 q = m / n;
        steps += q;
        i64 t = m - q * n;
        m = n;
        n = t;
    }
    return steps;
}

typedef struct { i128 c; i64 s; } entry;

static int cmp_entry(const void *a, const void *b) {
    const entry *ea = (const entry *)a;
    const entry *eb = (const entry *)b;
    if (ea->c < eb->c) return -1;
    if (ea->c > eb->c) return 1;
    return 0;
}

static i64 compute_F_for_powers(int k) {
    i64 a = 1, b = 1;
    for (int i = 0; i < k; i++) { a *= 6; b *= 10; }

    i64 p2[32], p3[32], p5[32];
    p2[0] = 1; p3[0] = 1; p5[0] = 1;
    for (int i = 0; i < k; i++) {
        p2[i + 1] = p2[i] * 2;
        p3[i + 1] = p3[i] * 3;
        p5[i + 1] = p5[i] * 5;
    }

    int na = (k + 1) * (k + 1);
    i64 div_a[400], ua[400];
    int ia = 0;
    for (int i = 0; i <= k; i++)
        for (int j = 0; j <= k; j++) {
            div_a[ia] = p2[i] * p3[j];
            ua[ia] = a / div_a[ia];
            ia++;
        }

    int nb = (k + 1) * (k + 1);
    i64 div_b[400], vb[400];
    int ib = 0;
    for (int i = 0; i <= k; i++)
        for (int j = 0; j <= k; j++) {
            div_b[ib] = p2[i] * p5[j];
            vb[ib] = b / div_b[ib];
            ib++;
        }

    int max_entries = na * nb * 2;
    entry *arr = (entry *)malloc(sizeof(entry) * max_entries);
    int cnt = 0;

    for (int i = 0; i < na; i++) {
        i64 x = div_a[i], u = ua[i];
        for (int j = 0; j < nb; j++) {
            i64 y = div_b[j], v = vb[j];
            i64 s = euclid_subtraction_steps(x, y);

            i128 c1 = (i128)(x + y) * (i128)(u + v);
            arr[cnt].c = c1; arr[cnt].s = s; cnt++;

            i128 c2 = (i128)(x - y) * (i128)(u - v);
            if (c2 > 0) {
                i64 s2 = s - 1;
                if (s2 > 0) {
                    arr[cnt].c = c2; arr[cnt].s = s2; cnt++;
                }
            }
        }
    }

    qsort(arr, cnt, sizeof(entry), cmp_entry);

    i64 sum = 0;
    i128 prev_c = 0;
    i64 best_s = 0;
    int have = 0;
    for (int i = 0; i < cnt; i++) {
        if (have && arr[i].c == prev_c) {
            if (arr[i].s < best_s) best_s = arr[i].s;
        } else {
            if (have) sum += best_s;
            prev_c = arr[i].c;
            best_s = arr[i].s;
            have = 1;
        }
    }
    if (have) sum += best_s;

    free(arr);
    return sum;
}

long long p876_native(void) {
    i64 total = 0;
    for (int k = 1; k <= 18; k++) {
        total += compute_F_for_powers(k);
    }
    return total;
}
