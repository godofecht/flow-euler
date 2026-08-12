/* Project Euler 924: Larger Digit Permutation II
 *
 * B(n) = smallest number larger than n formed by rearranging digits of n, or 0.
 * a_0 = 0, a_n = a_{n-1}^2 + 2.
 * U(N) = sum_{n=1..N} B(a_n).
 * Answer: U(10^16) mod 1_000_000_007.
 *
 * Algorithm:
 *   U(N) = sum a_n + sum (B(a_n) - a_n)  (mod MOD)
 *   - sum_a_mod: cycle detection for a_n mod MOD
 *   - delta_small: direct computation for n=1..5
 *   - delta10: 10-digit next-perm deltas for n>=6
 *   - delta_bad_11: 11-digit next-perm for "bad" indices
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;
typedef __int128 i128;

#define MOD 1000000007LL

/* ---- Hash table for cycle detection ---- */

#define HT_SIZE (1 << 17)  /* 131072, enough for ~61K entries */
#define HT_MASK (HT_SIZE - 1)

typedef struct { i64 key; i64 val; } HTEntry;

static HTEntry ht_table[HT_SIZE];

static void ht_init(void) {
    memset(ht_table, 0, sizeof(ht_table));
    /* Use 0 as sentinel for empty (key=0, val=0 is valid for a_0=0 at index 0) */
    /* We handle this by initializing a_0=0 at index 0 before any inserts */
}

static i64 ht_get(i64 key) {
    i64 h = (key * 2654435761LL) & HT_MASK;
    if (h < 0) h = -h;
    for (int i = 0; i < HT_SIZE; i++) {
        i64 idx = (h + i) & HT_MASK;
        if (ht_table[idx].key == key && ht_table[idx].val != 0) {
            return ht_table[idx].val;
        }
        if (ht_table[idx].val == 0 && ht_table[idx].key == 0) {
            /* empty slot (but key=0 at val=0 could be a_0) */
            /* This is ambiguous; handle a_0=0 separately */
            return -1; /* not found */
        }
    }
    return -1;
}

static void ht_put(i64 key, i64 val) {
    i64 h = (key * 2654435761LL) & HT_MASK;
    if (h < 0) h = -h;
    for (int i = 0; i < HT_SIZE; i++) {
        i64 idx = (h + i) & HT_MASK;
        if (ht_table[idx].val == 0) {
            ht_table[idx].key = key;
            ht_table[idx].val = val;
            return;
        }
    }
}

/* ---- Next permutation (for small exact integers) ---- */

static i64 next_perm_inplace(char *digs, int len) {
    int i = len - 2;
    while (i >= 0 && digs[i] >= digs[i + 1]) i--;
    if (i < 0) return 0; /* no next permutation */
    int j = len - 1;
    while (digs[j] <= digs[i]) j--;
    char tmp = digs[i]; digs[i] = digs[j]; digs[j] = tmp;
    int l = i + 1, r = len - 1;
    while (l < r) {
        tmp = digs[l]; digs[l] = digs[r]; digs[r] = tmp;
        l++; r--;
    }
    return 1;
}

static i64 B_small(i64 n) {
    char digs[20];
    int len = 0;
    i64 t = n;
    if (t == 0) { digs[0] = 0; len = 1; }
    else { while (t > 0) { digs[len++] = t % 10; t /= 10; }
           /* reverse */
           for (int i = 0; i < len/2; i++) { char tmp = digs[i]; digs[i] = digs[len-1-i]; digs[len-1-i] = tmp; } }
    if (!next_perm_inplace(digs, len)) return 0;
    i64 y = 0;
    for (int i = 0; i < len; i++) y = y * 10 + digs[i];
    return y;
}

/* ---- Fixed-width next permutation (leading zeros allowed) ---- */

static i64 next_perm_fixed(i64 x, int k, char *buf) {
    /* Extract k digits */
    i64 t = x;
    for (int i = k - 1; i >= 0; i--) {
        buf[i] = t % 10;
        t /= 10;
    }
    /* Find pivot */
    int i = k - 2;
    while (i >= 0 && buf[i] >= buf[i + 1]) i--;
    if (i < 0) return -1; /* no next permutation */
    int j = k - 1;
    while (buf[j] <= buf[i]) j--;
    char tmp = buf[i]; buf[i] = buf[j]; buf[j] = tmp;
    int l = i + 1, r = k - 1;
    while (l < r) {
        tmp = buf[l]; buf[l] = buf[r]; buf[r] = tmp;
        l++; r--;
    }
    i64 y = 0;
    for (int d = 0; d < k; d++) y = y * 10 + buf[d];
    return y;
}

/* ---- sum_a_mod: cycle detection ---- */

static i64 sum_a_mod(i64 N) {
    /* states[i] = a_i mod MOD, pref[i] = sum_{t=1..i} a_t mod MOD */
    i64 *states = (i64 *)malloc(200000 * sizeof(i64));
    i64 *pref = (i64 *)malloc(200000 * sizeof(i64));
    int n_states = 0;

    ht_init();
    states[0] = 0;
    ht_put(0, 0); /* a_0 = 0 at index 0 */
    n_states = 1;

    i64 x = 0;
    i64 mu = 0, lam = 0;

    while (1) {
        i64 nxt = (x * x + 2) % MOD;
        i64 idx = n_states;

        /* Check if nxt is in seen */
        i64 found = -1;
        if (nxt == 0) {
            found = 0; /* a_0 = 0 */
        } else {
            i64 h = (nxt * 2654435761LL) & HT_MASK;
            if (h < 0) h = -h;
            for (int i = 0; i < HT_SIZE; i++) {
                i64 pos = (h + i) & HT_MASK;
                if (ht_table[pos].val == 0) break; /* empty */
                if (ht_table[pos].key == nxt) { found = ht_table[pos].val; break; }
            }
        }

        if (found >= 0) {
            mu = found;
            lam = idx - found;
            break;
        }
        ht_put(nxt, idx);
        states[n_states++] = nxt;
        x = nxt;
    }

    /* Build prefix sums */
    pref[0] = 0;
    for (int i = 1; i < n_states; i++) {
        pref[i] = (pref[i - 1] + states[i]) % MOD;
    }

    if (N < n_states) {
        i64 result = pref[N];
        free(states); free(pref);
        return result;
    }

    i64 base_before = (mu > 0) ? pref[mu - 1] : 0;
    i64 cycle_sum = (pref[mu + lam - 1] - base_before) % MOD;
    if (cycle_sum < 0) cycle_sum += MOD;

    i64 cycle_terms = N - mu + 1;
    i64 full = cycle_terms / lam;
    i64 rem = cycle_terms % lam;

    i64 total = (base_before + (full % MOD) * cycle_sum) % MOD;
    if (rem) {
        i64 extra = (pref[mu + (int)rem - 1] - base_before) % MOD;
        if (extra < 0) extra += MOD;
        total = (total + extra) % MOD;
    }

    free(states); free(pref);
    return total;
}

/* ---- delta_small: n=1..5 ---- */

static i64 delta_small(i64 N) {
    i64 a = 0;
    i64 s = 0;
    int limit = (N < 5) ? (int)N : 5;
    for (int n = 1; n <= limit; n++) {
        a = a * a + 2;
        i64 b = B_small(a);
        i64 d = (b - a) % MOD;
        if (d < 0) d += MOD;
        s = (s + d) % MOD;
    }
    return s;
}

/* ---- delta10_and_bad ---- */

static i64 mod128(i128 x, i128 m) {
    i128 r = x % m;
    if (r < 0) r += m;
    return (i64)r;
}

static void delta10_and_bad(i64 N, i64 *out_total, i64 *out_first_bad_n, i64 *out_step) {
    *out_total = 0;
    *out_first_bad_n = -1;
    *out_step = 0;

    if (N <= 5) return;

    int k = 10;
    i128 m = 1;
    for (int i = 0; i < k; i++) m *= 10;
    i64 step = 8;
    for (int i = 0; i < k - 2; i++) step *= 5;
    *out_step = step;

    /* Compute a_6 mod 10^10 */
    i128 x = 0;
    for (int i = 0; i < 6; i++) x = (x * x + 2) % m;
    i128 start = x;

    i64 total_terms = N - 5;
    i64 q = total_terms / step;
    i64 r = total_terms % step;

    char buf[16];
    i64 cycle_sum = 0;
    i64 rem_sum = 0;
    i64 bad_step = -1;

    for (i64 i = 1; i <= step; i++) {
        i64 y = next_perm_fixed((i64)x, k, buf);
        if (y < 0) {
            if (bad_step != -1) {
                fprintf(stderr, "More than one bad position!\n");
                exit(1);
            }
            bad_step = i;
        } else {
            i64 d = y - (i64)x;
            cycle_sum = (cycle_sum + d) % MOD;
            if (cycle_sum < 0) cycle_sum += MOD;
            if (i <= r) {
                rem_sum = (rem_sum + d) % MOD;
                if (rem_sum < 0) rem_sum += MOD;
            }
        }
        x = (x * x + 2) % m;
    }

    if (x != start) {
        fprintf(stderr, "Cycle verification failed!\n");
        exit(1);
    }
    if (bad_step < 0) {
        fprintf(stderr, "No bad step found!\n");
        exit(1);
    }

    *out_first_bad_n = bad_step + 5;
    *out_total = ((q % MOD) * cycle_sum % MOD + rem_sum) % MOD;
    if (*out_total < 0) *out_total += MOD;
}

/* ---- delta_bad_11 ---- */

static i64 delta_bad_11(i64 N, i64 first_bad_n, i64 step) {
    if (first_bad_n < 0 || N < first_bad_n) return 0;

    i128 m11 = 1;
    for (int i = 0; i < 11; i++) m11 *= 10;

    i64 targets[5];
    for (int t = 0; t < 5; t++) targets[t] = first_bad_n + t * step;
    i64 max_n = targets[4];

    /* Simulate a_n mod 10^11 up to max_n */
    i128 x = 0;
    i64 vals[5];
    int idx = 0;
    for (i64 n = 1; n <= max_n; n++) {
        x = (x * x + 2) % m11;
        if (n == targets[idx]) {
            vals[idx] = (i64)x;
            idx++;
            if (idx == 5) break;
        }
    }

    char buf[16];
    i64 deltas[5];
    for (int i = 0; i < 5; i++) {
        i64 y = next_perm_fixed(vals[i], 11, buf);
        if (y < 0) {
            fprintf(stderr, "11-digit perm also None!\n");
            exit(1);
        }
        deltas[i] = (y - vals[i]) % MOD;
        if (deltas[i] < 0) deltas[i] += MOD;
    }

    /* Count bad indices <= N */
    i64 T = 1 + (N - first_bad_n) / step;
    i64 base = T / 5;
    i64 rem = T % 5;

    i64 total = 0;
    for (int i = 0; i < 5; i++) {
        i64 c = base + (i < rem ? 1 : 0);
        total = (total + (c % MOD) * deltas[i]) % MOD;
    }
    return total;
}

/* ---- solve ---- */

long long p924_native(void) {
    i64 N = 1;
    for (int i = 0; i < 16; i++) N *= 10;

    i64 s_a = sum_a_mod(N);
    i64 d_small = delta_small(N);

    i64 d10, first_bad_n, step;
    delta10_and_bad(N, &d10, &first_bad_n, &step);
    i64 d_bad = delta_bad_11(N, first_bad_n, step);

    i64 result = (s_a + d_small + d10 + d_bad) % MOD;
    if (result < 0) result += MOD;
    return result;
}
