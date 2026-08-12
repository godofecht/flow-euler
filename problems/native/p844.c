// Project Euler 844: k-Markov Numbers
// S(10^18, 10^18) mod 1405695061
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

#define MOD 1405695061LL

static i64 mmul(i64 a, i64 b) { return (i64)((i128)a * b % MOD); }
static i64 madd(i64 a, i64 b) { return (a + b) % MOD; }
static i64 msub(i64 a, i64 b) { return (a - b % MOD + MOD) % MOD; }

static i64 modpow(i64 base, i64 exp, i64 mod) {
    i64 r = 1; base %= mod;
    while (exp > 0) {
        if (exp & 1) r = (i64)((i128)r * base % mod);
        base = (i64)((i128)base * base % mod);
        exp >>= 1;
    }
    return r;
}

// ---- integer square root ----
static i64 isqrt64(i64 n) {
    if (n <= 0) return 0;
    i64 x = (i64)sqrtl((long double)n);
    while (x > 0 && x * x > n) x--;
    while ((x + 1) * (x + 1) <= n) x++;
    return x;
}

// ---- polynomial helpers for large-k regions ----
static i64 msum1(i64 n) {
    if (n <= 0) return 0;
    i64 a = n % MOD, b = (n + 1) % MOD;
    return mmul(mmul(a, b), modpow(2, MOD - 2, MOD));
}
static i64 msum2(i64 n) {
    if (n <= 0) return 0;
    i64 a = n % MOD, b = (n + 1) % MOD, c = (2 * n + 1) % MOD;
    return mmul(mmul(mmul(a, b), c), modpow(6, MOD - 2, MOD));
}
static i64 msum3(i64 n) {
    if (n <= 0) return 0;
    i64 s = msum1(n);
    return mmul(s, s);
}
static i64 mrange_sum(i64 (*func)(i64), i64 l, i64 r) {
    if (l > r) return 0;
    return msub(func(r), func(l - 1));
}

// ---- closed-form k helpers ----
static i128 poly_m3(i64 k) { return (i128)k*k*k*k - 2*(i128)k*k*k + k - 1; }
static i128 poly_a2(i64 k) { return (i128)k*k - k - 1; }
static i128 poly_a3(i64 k) { return (i128)k*k*k - k*k - 2*k + 1; }

static i64 max_k_monotone(i64 limit, i64 start, int (*pred)(i64)) {
    if (limit < start) return limit;
    if (!pred(start)) return start - 1;
    i64 lo = start, hi = start;
    while (hi < limit && pred(hi)) {
        lo = hi;
        hi = (hi * 2 > limit) ? limit : hi * 2;
    }
    if (pred(hi)) return hi;
    while (lo + 1 < hi) {
        i64 mid = (lo + hi) / 2;
        if (pred(mid)) lo = mid;
        else hi = mid;
    }
    return lo;
}

static i64 N_val;
static int pred_m3(i64 k) { return poly_m3(k) <= (i128)N_val; }
static int pred_a3(i64 k) { return poly_a3(k) <= (i128)N_val; }

static i64 max_k_three_nonones(i64 N, i64 K) {
    if (K < 3) return K;
    N_val = N;
    return max_k_monotone(K, 3, pred_m3);
}

static i64 max_k_a2(i64 N, i64 K) {
    if (K < 1) return 0;
    i64 disc = 1 + 4*(N+1);
    i64 r = (1 + isqrt64(disc)) / 2;
    if (r > K) r = K;
    while (r > 0 && poly_a2(r) > N) r--;
    return r;
}

static i64 max_k_a3(i64 N, i64 K) {
    if (K < 3) return K;
    if (poly_a3(3) > N) return 2;
    N_val = N;
    return max_k_monotone(K, 3, pred_a3);
}

// ---- Mk_sum via DFS ----
// State: sorted non-one values (up to 8), product as i128
#define MAX_NONONES 8
#define HT_CAP 2048

typedef struct {
    i64 vals[MAX_NONONES];
    int n;
} State;

typedef struct {
    State keys[HT_CAP];
    char used[HT_CAP];
} StateSet;

static u64 hash_state(const State *s) {
    u64 h = 1469598103934665603ULL;
    for (int i = 0; i < s->n; i++) {
        h ^= (u64)s->vals[i];
        h *= 1099511628211ULL;
    }
    return h;
}

static int state_eq(const State *a, const State *b) {
    if (a->n != b->n) return 0;
    for (int i = 0; i < a->n; i++)
        if (a->vals[i] != b->vals[i]) return 0;
    return 1;
}

static void state_clear(StateSet *set) {
    memset(set->used, 0, HT_CAP);
}

static int state_contains(StateSet *set, const State *s) {
    u64 h = hash_state(s) % HT_CAP;
    while (set->used[h]) {
        if (state_eq(&set->keys[h], s)) return 1;
        h = (h + 1) % HT_CAP;
    }
    return 0;
}

static void state_insert(StateSet *set, const State *s) {
    u64 h = hash_state(s) % HT_CAP;
    while (set->used[h]) {
        if (state_eq(&set->keys[h], s)) return;
        h = (h + 1) % HT_CAP;
    }
    set->used[h] = 1;
    set->keys[h] = *s;
}

// Seen-numbers hash set for i64 values
#define SN_CAP 4096
typedef struct {
    i64 keys[SN_CAP];
    char used[SN_CAP];
} NumSet;

static void num_clear(NumSet *set) {
    memset(set->used, 0, SN_CAP);
}

static int num_contains(NumSet *set, i64 v) {
    u64 h = (u64)v % SN_CAP;
    while (set->used[h]) {
        if (set->keys[h] == v) return 1;
        h = (h + 1) % SN_CAP;
    }
    return 0;
}

static void num_insert(NumSet *set, i64 v) {
    u64 h = (u64)v % SN_CAP;
    while (set->used[h]) {
        if (set->keys[h] == v) return;
        h = (h + 1) % SN_CAP;
    }
    set->used[h] = 1;
    set->keys[h] = v;
}

// Stack entry
typedef struct {
    State state;
    i128 prod;
} StackEntry;

static i64 mk_sum(i64 k, i64 N, i64 mod) {
    if (N < 1) return 0;
    i64 s = 1 % mod;
    if (k - 1 > N) return s;

    StateSet visited;
    state_clear(&visited);
    NumSet seen;
    num_clear(&seen);
    num_insert(&seen, 1);

    StackEntry stack[4096];
    int sp = 0;

    State start = { .n = 0 };
    stack[sp].state = start;
    stack[sp].prod = 1;
    sp++;
    state_insert(&visited, &start);

    while (sp > 0) {
        sp--;
        State non_ones = stack[sp].state;
        i128 prod_non_ones = stack[sp].prod;
        int ones = k - non_ones.n;

        // Record coordinates
        for (int i = 0; i < non_ones.n; i++) {
            i64 v = non_ones.vals[i];
            if (v <= N && !num_contains(&seen, v)) {
                num_insert(&seen, v);
                s = madd(s, v % mod);
            }
        }

        // Distinct values to jump
        // Collect unique values from non_ones plus 1 if ones > 0
        i64 jump_vals[MAX_NONONES + 1];
        int n_jump = 0;
        for (int i = 0; i < non_ones.n; i++) {
            if (i == 0 || non_ones.vals[i] != non_ones.vals[i-1]) {
                jump_vals[n_jump++] = non_ones.vals[i];
            }
        }
        if (ones > 0) {
            jump_vals[n_jump++] = 1;
        }

        for (int vi = 0; vi < n_jump; vi++) {
            i64 v = jump_vals[vi];
            i128 new_val;
            State new_state;
            i128 new_prod;

            if (v == 1) {
                if (ones == 0) continue;
                // Check overflow: if prod > (N+1)/k, skip
                if (prod_non_ones > (i128)(N + 1) / k) continue;
                new_val = (i128)k * prod_non_ones - 1;
                if (new_val <= 1 || new_val > N) continue;
                // Build new state: insert new_val into sorted non_ones
                new_state = non_ones;
                int pos = new_state.n;
                for (int j = 0; j < new_state.n; j++) {
                    if (new_val < new_state.vals[j]) {
                        pos = j;
                        break;
                    }
                }
                for (int j = new_state.n; j > pos; j--)
                    new_state.vals[j] = new_state.vals[j-1];
                new_state.vals[pos] = (i64)new_val;
                new_state.n++;
                new_prod = prod_non_ones * new_val;
            } else {
                // Jump a non-one: new = k*(prod/v) - v
                i128 prod_div_v = prod_non_ones / v;
                // Check overflow
                if (prod_div_v > (i128)(N + v) / k) continue;
                new_val = (i128)k * prod_div_v - v;
                if (new_val <= v || new_val > N) continue;
                // Replace v with new_val in sorted order
                new_state = non_ones;
                // Find and remove v, insert new_val
                int found = -1;
                for (int j = 0; j < new_state.n; j++) {
                    if (new_state.vals[j] == v) { found = j; break; }
                }
                if (found < 0) continue; // shouldn't happen
                // Remove v
                for (int j = found; j < new_state.n - 1; j++)
                    new_state.vals[j] = new_state.vals[j+1];
                new_state.n--;
                // Insert new_val
                int pos = new_state.n;
                for (int j = 0; j < new_state.n; j++) {
                    if (new_val < new_state.vals[j]) {
                        pos = j;
                        break;
                    }
                }
                for (int j = new_state.n; j > pos; j--)
                    new_state.vals[j] = new_state.vals[j-1];
                new_state.vals[pos] = (i64)new_val;
                new_state.n++;
                new_prod = prod_div_v * new_val;
            }

            if (!state_contains(&visited, &new_state)) {
                state_insert(&visited, &new_state);
                stack[sp].state = new_state;
                stack[sp].prod = new_prod;
                sp++;
            }
        }
    }

    return s % mod;
}

long long p844_native(void) {
    i64 K = 1000000000000000000LL; // 10^18
    i64 N = 1000000000000000000LL;
    i64 mod = MOD;

    if (K < 3 || N < 1) return 0;

    i64 K_eff = (K < N + 1) ? K : N + 1;
    i64 total = 0;

    i64 cutoff = (K_eff < max_k_three_nonones(N, K_eff)) ? K_eff : max_k_three_nonones(N, K_eff);

    // Enumerate exactly for k <= cutoff
    for (i64 k = 3; k <= cutoff; k++) {
        total = madd(total, mk_sum(k, N, mod));
    }

    i64 start = cutoff + 1;
    if (start > K_eff) {
        if (K > K_eff) total = madd(total, (K - K_eff) % mod);
        return total % mod;
    }

    i64 k3 = (K_eff < max_k_a3(N, K_eff)) ? K_eff : max_k_a3(N, K_eff);
    i64 k2 = (K_eff < max_k_a2(N, K_eff)) ? K_eff : max_k_a2(N, K_eff);

    // Region 1: a3 <= N => M_k = k^3 - 2k
    {
        i64 l = start, r = k3;
        if (l <= r) {
            i64 part = msub(mrange_sum(msum3, l, r), mmul(2, mrange_sum(msum1, l, r)));
            total = madd(total, part);
        }
    }

    // Region 2: a2 <= N < a3 => M_k = k^2 - 1
    {
        i64 l = (start > k3 + 1) ? start : k3 + 1;
        i64 r = k2;
        if (l <= r) {
            i64 cnt = (r - l + 1) % mod;
            i64 part = msub(mrange_sum(msum2, l, r), cnt);
            total = madd(total, part);
        }
    }

    // Region 3: a2 > N => M_k = k
    {
        i64 l = (start > k2 + 1) ? start : k2 + 1;
        i64 r = K_eff;
        if (l <= r) {
            i64 part = mrange_sum(msum1, l, r);
            total = madd(total, part);
        }
    }

    // k > N+1 contributes 1 each
    if (K > K_eff) {
        total = madd(total, (K - K_eff) % mod);
    }

    return total % mod;
}
