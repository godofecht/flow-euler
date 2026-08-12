/* Project Euler 943 - Kolakoski-like sequence counts.
   Port of the Python reference solver. Uses a memoized recursive
   solver with an open-addressing hash table (generation-stamped)
   to count a-runs and b-runs in the generalized Kolakoski sequence. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef struct {
    uint64_t count_a;
    uint64_t count_b;
    uint64_t next_state;
} Result;

typedef struct {
    uint64_t key;       /* 0 means empty (key 0 never used because cache keys are large) */
    uint64_t gen;       /* generation stamp */
    Result    res;
} Entry;

typedef struct {
    Entry *table;
    uint64_t mask;
    uint64_t gen;       /* current generation */
    uint64_t count;     /* number of entries in current generation */
    int a, b;
} Solver;

static Solver g_solver;

static void solver_init(int a, int b) {
    g_solver.a = a;
    g_solver.b = b;
    g_solver.gen++;
    g_solver.count = 0;
    /* table starts at 2^16; grows if needed */
}

static void solver_ensure_table(uint64_t need) {
    if (g_solver.table == NULL) {
        uint64_t sz = 1ULL << 16;
        while (sz < need * 3) sz <<= 1;  /* load factor < 0.33 */
        g_solver.table = calloc(sz, sizeof(Entry));
        g_solver.mask = sz - 1;
    } else {
        uint64_t cap = g_solver.mask + 1;
        if (g_solver.count * 3 >= cap) {
            /* resize */
            uint64_t newcap = cap;
            while (newcap < g_solver.count * 6) newcap <<= 1;
            Entry *old = g_solver.table;
            uint64_t oldcap = cap;
            g_solver.table = calloc(newcap, sizeof(Entry));
            g_solver.mask = newcap - 1;
            for (uint64_t i = 0; i < oldcap; i++) {
                if (old[i].gen == g_solver.gen && old[i].key != 0) {
                    uint64_t h = (old[i].key * 2654435761ULL) & g_solver.mask;
                    while (g_solver.table[h].key != 0) h = (h + 1) & g_solver.mask;
                    g_solver.table[h] = old[i];
                }
            }
            free(old);
        }
    }
}

static Entry *solver_lookup(uint64_t key) {
    uint64_t h = (key * 2654435761ULL) & g_solver.mask;
    for (;;) {
        Entry *e = &g_solver.table[h];
        if (e->key == 0 || e->gen != g_solver.gen) return e;  /* empty slot */
        if (e->key == key) return e;  /* found */
        h = (h + 1) & g_solver.mask;
    }
}

static Result calc(uint64_t state, int level, uint64_t maxlen) {
    if (maxlen == 0) {
        Result r = {0, 0, state};
        return r;
    }

    uint64_t length_bit = (uint64_t)2 << level;
    uint64_t bit = state & length_bit;
    int run_len = bit ? g_solver.b : g_solver.a;
    uint64_t count = (uint64_t)run_len < maxlen ? (uint64_t)run_len : maxlen;

    if (level == 0) {
        if ((state & 1) == 0) {
            Result r = {count, 0, state ^ 1};
            return r;
        } else {
            Result r = {0, count, state ^ 1};
            return r;
        }
    }

    uint64_t produced_a = 0, produced_b = 0;
    uint64_t substate = state ^ bit;

    for (uint64_t i = 0; i < count; i++) {
        uint64_t child_key = substate + ((uint64_t)2 << level);
        solver_ensure_table(g_solver.count + 1);
        Entry *e = solver_lookup(child_key);
        Result child;
        if (e->key == child_key && e->gen == g_solver.gen) {
            child = e->res;
            uint64_t child_total = child.count_a + child.count_b;
            if (produced_a + produced_b + child_total > maxlen) {
                child = calc(substate, level - 1, maxlen - produced_a - produced_b);
            }
        } else {
            child = calc(substate, level - 1, maxlen - produced_a - produced_b);
        }
        produced_a += child.count_a;
        produced_b += child.count_b;
        substate = child.next_state;
    }

    Result res;
    res.count_a = produced_a;
    res.count_b = produced_b;
    res.next_state = substate ^ bit ^ ((uint64_t)1 << level);

    /* store in cache */
    uint64_t cache_key = state + ((uint64_t)4 << level);
    solver_ensure_table(g_solver.count + 1);
    Entry *e = solver_lookup(cache_key);
    if (e->key != cache_key || e->gen != g_solver.gen) {
        e->key = cache_key;
        e->gen = g_solver.gen;
        g_solver.count++;
    }
    e->res = res;
    return res;
}

static uint64_t compute_T(int a, int b, uint64_t limit) {
    solver_init(a, b);
    int level = 0;
    Result res = {0, 0, 0};
    for (;;) {
        res = calc(0, level, limit);
        level++;
        if (res.count_a + res.count_b >= limit || level >= 64) break;
    }
    return res.count_a * (uint64_t)a + res.count_b * (uint64_t)b;
}

long long p943_native(void) {
    /* self-tests */
    if (compute_T(2, 3, 10) != 25) {
        fprintf(stderr, "self-test 1 failed\n");
        return -1;
    }
    if (compute_T(4, 2, 10000) != 30004) {
        fprintf(stderr, "self-test 2 failed\n");
        return -1;
    }

    uint64_t MOD = 2233222333ULL;
    uint64_t N = 22332223332233ULL;
    uint64_t total = 0;
    for (int a = 2; a < 224; a++) {
        for (int b = 2; b < 224; b++) {
            if (a == b) continue;
            uint64_t contribution = compute_T(a, b, N) % MOD;
            total = (total + contribution) % MOD;
        }
    }
    return (long long)(total % MOD);
}
