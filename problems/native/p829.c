// Project Euler 829: Integral Fusion
// sum_{n=2..31} M(n) where M(n) is smallest integer whose factor tree
// has the same shape as T(n!!).
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef unsigned long long u64;
typedef __int128 i128;

// ---- Miller-Rabin ----
static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((i128)a * b % m); }
static u64 powmod(u64 a, u64 d, u64 n) {
    u64 r = 1; a %= n;
    while (d) { if (d&1) r = mulmod(r, a, n); a = mulmod(a, a, n); d >>= 1; }
    return r;
}
static int is_prime_u64(u64 n) {
    if (n < 2) return 0;
    u64 sp[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++) { if (n == sp[i]) return 1; if (n % sp[i] == 0) return 0; }
    u64 d = n-1; int s = 0;
    while (d%2==0) { d/=2; s++; }
    u64 bases[] = {2,325,9375,28178,450775,9780504,1795265022};
    for (int i = 0; i < 7; i++) {
        u64 a = bases[i] % n; if (a == 0) continue;
        u64 x = powmod(a, d, n);
        if (x == 1 || x == n-1) continue;
        int found = 0;
        for (int j = 0; j < s-1; j++) { x = mulmod(x, x, n); if (x == n-1) { found = 1; break; } }
        if (!found) return 0;
    }
    return 1;
}

static u64 gcd_u64(u64 a, u64 b) { while (b) { u64 t = a%b; a=b; b=t; } return a; }

static u64 pollard_rho_u64(u64 n) {
    if (n % 2 == 0) return 2;
    if (n % 3 == 0) return 3;
    u64 c_vals[] = {1,3,5,7,11,13,17,19,23};
    for (int ci = 0; ci < 9; ci++) {
        u64 c = c_vals[ci], y = 2, x, ys, g = 1, q = 1, r = 1, m = 128;
        while (g == 1) {
            x = y;
            for (u64 i = 0; i < r; i++) y = (mulmod(y, y, n) + c) % n;
            u64 k = 0;
            while (k < r && g == 1) {
                ys = y;
                u64 lim = r - k; if (m < lim) lim = m;
                for (u64 i = 0; i < lim; i++) {
                    y = (mulmod(y, y, n) + c) % n;
                    u64 diff = x > y ? x - y : y - x;
                    q = mulmod(q, diff, n);
                }
                g = gcd_u64(q, n);
                k += m;
            }
            r <<= 1;
        }
        if (g == n) { g = 1; y = ys; while (g == 1) { y = (mulmod(y, y, n) + c) % n; u64 diff = x>y?x-y:y-x; g = gcd_u64(diff, n); } }
        if (1 < g && g < n) return g;
    }
    for (u64 d = 5; d <= 1000000 && d*d <= n; d += 2) if (n % d == 0) return d;
    return n;
}

typedef struct { u64 p; int e; } Factor;
static int factorize_u64(u64 n, Factor *out) {
    int nf = 0;
    u64 sp[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++) {
        if (n % sp[i] == 0) {
            int e = 0; while (n % sp[i] == 0) { n /= sp[i]; e++; }
            out[nf].p = sp[i]; out[nf].e = e; nf++;
        }
    }
    u64 stack[64]; int sp2 = 0;
    if (n > 1) stack[sp2++] = n;
    while (sp2 > 0) {
        u64 m = stack[--sp2];
        if (m == 1) continue;
        if (is_prime_u64(m)) {
            int found = 0;
            for (int i = 0; i < nf; i++) if (out[i].p == m) { out[i].e++; found = 1; break; }
            if (!found) { out[nf].p = m; out[nf].e = 1; nf++; }
            continue;
        }
        u64 d = pollard_rho_u64(m);
        if (d == m) {
            int found = 0;
            for (int i = 0; i < nf; i++) if (out[i].p == m) { out[i].e++; found = 1; break; }
            if (!found) { out[nf].p = m; out[nf].e = 1; nf++; }
            continue;
        }
        stack[sp2++] = d; stack[sp2++] = m / d;
    }
    return nf;
}

static u64 isqrt_u64(u64 n) {
    if (n == 0) return 0;
    u64 x = (u64)sqrt((double)n);
    while (x*x > n) x--; while ((x+1)*(x+1) <= n) x++; return x;
}

static int gen_divisors(Factor *fac, int nf, u64 *divs, int max_divs) {
    int nd = 1; divs[0] = 1;
    for (int i = 0; i < nf; i++) {
        int cur = nd; u64 pe = 1;
        for (int e = 1; e <= fac[i].e; e++) {
            pe *= fac[i].p;
            for (int j = 0; j < cur; j++) {
                if (nd >= max_divs) return nd;
                divs[nd++] = divs[j] * pe;
            }
        }
    }
    for (int i = 1; i < nd; i++) { u64 key = divs[i]; int j = i-1; while (j>=0 && divs[j]>key) { divs[j+1]=divs[j]; j--; } divs[j+1]=key; }
    return nd;
}

static u64 best_divisor_le_sqrt(u64 n) {
    Factor fac[20];
    int nf = factorize_u64(n, fac);
    u64 divs[20000];
    int nd = gen_divisors(fac, nf, divs, 20000);
    u64 root = isqrt_u64(n);
    u64 best = 1;
    for (int i = 0; i < nd; i++) if (divs[i] <= root && divs[i] > best) best = divs[i];
    return best;
}

// ---- Shape tree ----
typedef struct ShapeNode {
    int is_leaf;
    struct ShapeNode *left, *right;
} ShapeNode;

static ShapeNode *leaf_node = NULL;
static ShapeNode *new_shape(int is_leaf, ShapeNode *l, ShapeNode *r) {
    ShapeNode *n = malloc(sizeof(ShapeNode));
    n->is_leaf = is_leaf; n->left = l; n->right = r;
    return n;
}

static int shape_eq(ShapeNode *a, ShapeNode *b) {
    if (a == b) return 1;
    if (a->is_leaf != b->is_leaf) return 0;
    if (a->is_leaf) return 1;
    return shape_eq(a->left, b->left) && shape_eq(a->right, b->right);
}

#define SCACHE_SIZE 8192
typedef struct { u64 key; ShapeNode *val; } SCacheEntry;
static SCacheEntry scache[SCACHE_SIZE];
static int scache_count = 0;

static ShapeNode *get_shape(u64 n) {
    for (int i = 0; i < scache_count; i++) if (scache[i].key == n) return scache[i].val;
    ShapeNode *result;
    if (is_prime_u64(n)) {
        result = leaf_node;
    } else {
        u64 d = best_divisor_le_sqrt(n);
        u64 a = d, b = n / d;
        if (a > b) { u64 t = a; a = b; b = t; }
        result = new_shape(0, get_shape(a), get_shape(b));
    }
    if (scache_count < SCACHE_SIZE) { scache[scache_count].key = n; scache[scache_count].val = result; scache_count++; }
    return result;
}

static int count_leaves(ShapeNode *sh) {
    if (sh->is_leaf) return 1;
    return count_leaves(sh->left) + count_leaves(sh->right);
}

// Check if n's factor tree shape matches sh
static int shape_matches(u64 n, ShapeNode *sh) {
    if (sh->is_leaf) return is_prime_u64(n);
    if (is_prime_u64(n)) return 0;
    u64 d = best_divisor_le_sqrt(n);
    u64 a = d, b = n / d;
    if (a > b) { u64 t = a; a = b; b = t; }
    return shape_eq(get_shape(a), sh->left) && shape_eq(get_shape(b), sh->right);
}

static u64 next_prime(u64 p) {
    if (p < 2) return 2;
    if (p == 2) return 3;
    u64 c = p + 2;
    while (!is_prime_u64(c)) c += 2;
    return c;
}

// ---- Min-heap for (value, ...) pairs ----
typedef struct { u64 val; int data; } HeapEntry;

static void hpush(HeapEntry **h, int *size, int *cap, u64 val, int data) {
    if (*size >= *cap) {
        *cap *= 2;
        *h = realloc(*h, *cap * sizeof(HeapEntry));
    }
    int pos = (*size)++;
    (*h)[pos].val = val; (*h)[pos].data = data;
    while (pos > 0) {
        int p = (pos-1)/2;
        if ((*h)[p].val <= (*h)[pos].val) break;
        HeapEntry t = (*h)[p]; (*h)[p] = (*h)[pos]; (*h)[pos] = t;
        pos = p;
    }
}
static void hpop(HeapEntry *h, int *size, u64 *val, int *data) {
    *val = h[0].val; *data = h[0].data;
    (*size)--;
    if (*size > 0) {
        h[0] = h[*size];
        int pos = 0;
        while (1) {
            int l = 2*pos+1, r = 2*pos+2, s = pos;
            if (l < *size && h[l].val < h[s].val) s = l;
            if (r < *size && h[r].val < h[s].val) s = r;
            if (s == pos) break;
            HeapEntry t = h[pos]; h[pos] = h[s]; h[s] = t;
            pos = s;
        }
    }
}

// ---- Shape sequence: generate numbers with given shape in increasing order ----
// For leaf: primes in order
// For node(L, R): products a*b where shape(a)=L, shape(b)=R, a<=b, shape(a*b)=(L,R)
//   Use a min-heap over (a_index, b_index) pairs.

typedef struct {
    ShapeNode *sh;
    u64 maxval;
    // For leaf
    u64 next_prime_val;
    // For node
    int left_idx, right_idx;  // next indices to try
    HeapEntry *heap;
    int heap_size, heap_cap;
    // In-heap set (i,j) -> use a hash or simple array
    int *seen_pairs;  // encoded as i*100000+j
    int seen_count, seen_cap;
    // Generated values
    u64 *values;
    int val_count, val_cap;
} Seq;

static Seq **seqs = NULL;
static int seqs_count = 0;

static Seq *find_seq(ShapeNode *sh) {
    for (int i = 0; i < seqs_count; i++) if (seqs[i]->sh == sh) return seqs[i];
    return NULL;
}

static u64 get_value(ShapeNode *sh, int idx);

// Ensure right[j] >= x, return smallest j with right[j] >= x, or -1
static int ensure_right_ge(Seq *s, u64 x) {
    Seq *rs = find_seq(s->sh->right);
    if (!rs) return -1;
    if (rs->val_count == 0) {
        u64 v = get_value(s->sh->right, 0);
        if (v == 0 || v > s->maxval) return -1;
    }
    while (rs->values[rs->val_count - 1] < x) {
        u64 v = get_value(s->sh->right, rs->val_count);
        if (v == 0 || v > s->maxval) break;
    }
    if (rs->values[rs->val_count - 1] < x) return -1;
    int lo = 0, hi = rs->val_count - 1;
    while (lo < hi) { int mid = (lo+hi)/2; if (rs->values[mid] < x) lo = mid+1; else hi = mid; }
    return lo;
}

static int seen_contains(Seq *s, int i, int j) {
    for (int k = 0; k < s->seen_count; k++) {
        if (s->seen_pairs[k*2] == i && s->seen_pairs[k*2+1] == j) return 1;
    }
    return 0;
}
static void seen_add(Seq *s, int i, int j) {
    if (s->seen_count >= s->seen_cap) {
        s->seen_cap *= 2;
        s->seen_pairs = realloc(s->seen_pairs, s->seen_cap * 2 * sizeof(int));
    }
    s->seen_pairs[s->seen_count*2] = i;
    s->seen_pairs[s->seen_count*2+1] = j;
    s->seen_count++;
}

static void try_push(Seq *s, int i, int j) {
    if (seen_contains(s, i, j)) return;
    u64 a = get_value(s->sh->left, i);
    u64 b = get_value(s->sh->right, j);
    if (a == 0 || b == 0) return;
    if (a > b) return;
    u64 prod = a * b;
    if (prod > s->maxval) return;
    seen_add(s, i, j);
    hpush(&s->heap, &s->heap_size, &s->heap_cap, prod, i * 1000000 + j);
}

static u64 next_candidate(Seq *s) {
    // Start
    if (s->left_idx == 0 && s->heap_size == 0) {
        u64 a0 = get_value(s->sh->left, 0);
        if (a0 == 0) return 0;
        int j0 = ensure_right_ge(s, a0);
        if (j0 >= 0) try_push(s, 0, j0);
        s->left_idx = 1;
    }

    while (1) {
        if (s->heap_size == 0) {
            u64 a = get_value(s->sh->left, s->left_idx);
            if (a == 0) return 0;
            // Lower bound: a*a. If heap is empty and a*a > maxval, stop.
            if (a * a > s->maxval) {
                // But maybe a smaller b exists... no, b >= a, so a*b >= a*a
                return 0;
            }
            int j0 = ensure_right_ge(s, a);
            if (j0 < 0) return 0;
            try_push(s, s->left_idx, j0);
            s->left_idx++;
            continue;
        }
        u64 prod;
        int ij;
        hpop(s->heap, &s->heap_size, &prod, &ij);
        int i = ij / 1000000, j = ij % 1000000;

        // Maybe add more i values
        while (1) {
            u64 a = get_value(s->sh->left, s->left_idx);
            if (a == 0) break;
            if (a * a > prod && s->heap_size > 0) break;
            int j0 = ensure_right_ge(s, a);
            if (j0 >= 0) try_push(s, s->left_idx, j0);
            s->left_idx++;
        }

        // Advance j for this i
        try_push(s, i, j + 1);
        return prod;
    }
}

static u64 next_value(Seq *s) {
    if (s->sh->is_leaf) {
        u64 p = s->next_prime_val;
        if (p > s->maxval) return 0;
        s->next_prime_val = next_prime(p);
        return p;
    }
    u64 last = s->val_count > 0 ? s->values[s->val_count - 1] : 0;
    while (1) {
        u64 cand = next_candidate(s);
        if (cand == 0) return 0;
        if (last != 0 && cand == last) continue;
        if (shape_matches(cand, s->sh)) return cand;
    }
}

static u64 get_value(ShapeNode *sh, int idx) {
    Seq *s = find_seq(sh);
    if (!s) return 0;
    while (s->val_count <= idx) {
        u64 v = next_value(s);
        if (v == 0) return 0;
        if (s->val_count >= s->val_cap) {
            s->val_cap *= 2;
            s->values = realloc(s->values, s->val_cap * sizeof(u64));
        }
        s->values[s->val_count++] = v;
    }
    return s->values[idx];
}

// Collect all unique subshapes
static ShapeNode **all_shapes_arr = NULL;
static int all_shapes_arr_count = 0, all_shapes_arr_cap = 0;

static void collect_shape(ShapeNode *sh) {
    for (int i = 0; i < all_shapes_arr_count; i++) if (all_shapes_arr[i] == sh) return;
    if (all_shapes_arr_count >= all_shapes_arr_cap) {
        all_shapes_arr_cap = all_shapes_arr_cap ? all_shapes_arr_cap * 2 : 64;
        all_shapes_arr = realloc(all_shapes_arr, all_shapes_arr_cap * sizeof(ShapeNode *));
    }
    all_shapes_arr[all_shapes_arr_count++] = sh;
    if (!sh->is_leaf) { collect_shape(sh->left); collect_shape(sh->right); }
}

static u64 double_factorial(int n) {
    u64 r = 1;
    for (int k = n; k > 1; k -= 2) r *= k;
    return r;
}

long long p829_native(void) {
    leaf_node = new_shape(1, NULL, NULL);

    int max_n = 31;
    u64 maxval = double_factorial(max_n);

    ShapeNode *targets[32];
    for (int n = 2; n <= max_n; n++) targets[n] = get_shape(double_factorial(n));

    for (int n = 2; n <= max_n; n++) collect_shape(targets[n]);

    // Sort by leaf count
    for (int i = 1; i < all_shapes_arr_count; i++) {
        ShapeNode *key = all_shapes_arr[i];
        int j = i - 1;
        while (j >= 0 && count_leaves(all_shapes_arr[j]) > count_leaves(key)) {
            all_shapes_arr[j + 1] = all_shapes_arr[j]; j--;
        }
        all_shapes_arr[j + 1] = key;
    }

    // Create sequences
    seqs = malloc(all_shapes_arr_count * sizeof(Seq *));
    for (int i = 0; i < all_shapes_arr_count; i++) {
        Seq *s = calloc(1, sizeof(Seq));
        s->sh = all_shapes_arr[i];
        s->maxval = maxval;
        s->next_prime_val = 2;
        s->left_idx = 0;
        s->heap_cap = 4096;
        s->heap = malloc(s->heap_cap * sizeof(HeapEntry));
        s->heap_size = 0;
        s->seen_cap = 256;
        s->seen_pairs = malloc(s->seen_cap * 2 * sizeof(int));
        s->seen_count = 0;
        s->val_cap = 256;
        s->values = malloc(s->val_cap * sizeof(u64));
        s->val_count = 0;
        seqs[seqs_count++] = s;
    }

    u64 total = 0;
    for (int n = 2; n <= max_n; n++) {
        u64 m = get_value(targets[n], 0);
        total += m;
    }
    return (long long)total;
}
