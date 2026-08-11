/*
 * Project Euler 774 - Divide and Rule
 *
 * Counts sequences a_1..a_l with 0 <= a_i <= m such that every pair of
 * equal values a_i = a_j (i < j) has gcd(a_i, a_{i+1}, ..., a_j) > 1.
 *
 * Ported from the reference Python solver (solvers/774.py). The recursive
 * D(length, bound, left, right) is memoised with an open-addressing hash
 * table keyed on the four parameters.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MOD 998244353LL
#define TOP  (-1)
#define ODD  (-2)
#define MAX_N 123

/* ---- Fibonacci table (mod MOD) ---- */

static long long FIB[MAX_N + 2];

static void init_fib(void) {
    FIB[0] = 0;
    FIB[1] = 1 % MOD;
    for (int i = 2; i < MAX_N + 2; i++)
        FIB[i] = (FIB[i - 1] + FIB[i - 2]) % MOD;
}

static long long fib(int index) {
    if (index >= 0) return FIB[index];
    int k = -index;
    if (k & 1) return FIB[k];
    long long v = -FIB[k] % MOD;
    if (v < 0) v += MOD;
    return v;
}

/* ---- State helpers ---- */

static int boundary_is_impossible(int state) { return state == 0; }

static int satisfies(int state, int value) {
    if (state == TOP) return 1;
    if (state == ODD) return (value & 1) == 1;
    return (value & state) != 0;
}

static int phi_even(int state) {
    if (state == TOP) return TOP;
    if (state == ODD) return 0;
    return state / 2;
}

static int phi_odd(int state) {
    if (state == TOP || state == ODD) return TOP;
    if (state & 1) return TOP;
    return state / 2;
}

/* ---- Memoisation hash table ---- */
/* Key: (length, bound, left, right).  Value: result mod MOD. */

typedef struct {
    long long key;       /* packed */
    long long val;
    int used;
} Entry;

typedef struct {
    Entry *entries;
    size_t cap;
    size_t count;
} HashMap;

static HashMap memo;

static long long pack_key(int length, int bound, int left, int right) {
    /* length: 0..123 (7 bits), bound: up to ~1.3e8 (28 bits),
       left/right: signed ints (32 bits each).  Use 64-bit hash via
       combining; store the four raw ints in a 128-bit-ish key by
       packing into a single 64-bit value is too narrow, so instead
       we store the four ints in a parallel array indexed by slot.
       Simpler: keep a struct of 4 ints as the key. */
    /* Not used; we compare fields directly. */
    (void)length; (void)bound; (void)left; (void)right;
    return 0;
}

typedef struct {
    int length, bound, left, right;
} Key;

static Key *memo_keys;   /* parallel to memo.entries */

static unsigned long hash_key(const Key *k) {
    unsigned long h = 1469598103934665603UL;   /* FNV offset */
    h ^= (unsigned)k->length;  h *= 1099511628211UL;
    h ^= (unsigned)k->bound;   h *= 1099511628211UL;
    h ^= (unsigned)k->left;    h *= 1099511628211UL;
    h ^= (unsigned)k->right;   h *= 1099511628211UL;
    return h;
}

static int key_eq(const Key *a, const Key *b) {
    return a->length == b->length && a->bound == b->bound &&
           a->left == b->left && a->right == b->right;
}

static void memo_init(size_t cap) {
    memo.cap = cap;
    memo.count = 0;
    memo.entries = calloc(cap, sizeof(Entry));
    memo_keys = calloc(cap, sizeof(Key));
}

static void memo_insert_raw(Key k, long long v);

static void memo_grow(void) {
    Entry *old_e = memo.entries;
    Key *old_k = memo_keys;
    size_t old_cap = memo.cap;
    memo.cap *= 2;
    memo.entries = calloc(memo.cap, sizeof(Entry));
    memo_keys = calloc(memo.cap, sizeof(Key));
    memo.count = 0;
    for (size_t i = 0; i < old_cap; i++) {
        if (old_e[i].used) {
            Key k = old_k[i];
            long long v = old_e[i].val;
            memo_insert_raw(k, v);
        }
    }
    free(old_e);
    free(old_k);
}

static void memo_insert_raw(Key k, long long v) {
    unsigned long h = hash_key(&k) % memo.cap;
    while (memo.entries[h].used) {
        h = (h + 1) % memo.cap;
    }
    memo.entries[h].used = 1;
    memo.entries[h].val = v;
    memo_keys[h] = k;
    memo.count++;
}

/* Returns 1 if found (writes *out), 0 otherwise. */
static int memo_lookup(Key k, long long *out) {
    if (memo.cap == 0) return 0;
    unsigned long h = hash_key(&k) % memo.cap;
    while (memo.entries[h].used) {
        if (key_eq(&memo_keys[h], &k)) {
            *out = memo.entries[h].val;
            return 1;
        }
        h = (h + 1) % memo.cap;
    }
    return 0;
}

static void memo_put(Key k, long long v) {
    if (memo.count * 10 >= memo.cap * 7)
        memo_grow();
    memo_insert_raw(k, v);
}

/* ---- Core recursion ---- */

static long long D(int length, int bound, int left, int right);

static long long D(int length, int bound, int left, int right) {
    if (boundary_is_impossible(left) || boundary_is_impossible(right))
        return 0;

    if (length == 0)
        return (left == TOP && right == TOP) ? 1 : 0;

    if (bound <= 1) {
        if (length == 1) {
            long long s = 0;
            for (int value = 0; value <= bound; value++)
                if (satisfies(left, value) && satisfies(right, value))
                    s++;
            return s;
        }
        if (bound == 0)
            return 0;
        return (satisfies(left, 1) && satisfies(right, 1)) ? 1 : 0;
    }

    Key k = { length, bound, left, right };
    long long cached;
    if (memo_lookup(k, &cached)) return cached;

    long long total;

    if (bound % 2 == 0) {
        int marked_bound = bound;
        total = D(length, bound - 1, left, right);

        for (int split = 1; split <= length; split++) {
            int prefix_len = split - 1;
            int suffix_len = length - split;

            long long prefix_count;
            if (prefix_len == 0)
                prefix_count = satisfies(left, bound) ? 1 : 0;
            else
                prefix_count = D(prefix_len, bound - 1, left, marked_bound);
            if (prefix_count == 0)
                continue;

            long long suffix_count;
            if (suffix_len == 0)
                suffix_count = satisfies(right, bound) ? 1 : 0;
            else
                suffix_count = D(suffix_len, bound, marked_bound, right);

            total = (total + prefix_count * suffix_count) % MOD;
        }

        memo_put(k, total);
        return total;
    }

    /* bound is odd */
    int reduced_bound = (bound - 1) / 2;
    int left_even = phi_even(left);
    int left_odd = phi_odd(left);
    int right_even = phi_even(right);
    int right_odd = phi_odd(right);

    total = D(length, reduced_bound, left_even, right_even) * fib(length) % MOD;
    total = (total + D(length, reduced_bound, left_even, right_odd) * fib(length - 1)) % MOD;
    total = (total + D(length, reduced_bound, left_odd, right_even) * fib(length - 1)) % MOD;
    total = (total + D(length, reduced_bound, left_odd, right_odd) * fib(length - 2)) % MOD;
    total %= MOD;

    for (int cut = 1; cut < length; cut++) {
        long long prefix = D(cut, reduced_bound, left_odd, TOP) * fib(cut - 2) % MOD;
        if (cut > 1)
            prefix = (prefix + D(cut, reduced_bound, left_even, TOP) * fib(cut - 1)) % MOD;
        if (prefix) {
            total = (total + D(length - cut, bound, ODD, right) * prefix) % MOD;
        }
    }

    memo_put(k, total);
    return total;
}

static long long c_func(int length, int bound) {
    return D(length, bound, TOP, TOP);
}

long long p774_native(void) {
    init_fib();
    memo_init(1 << 16);
    return c_func(123, 123456789);
}
