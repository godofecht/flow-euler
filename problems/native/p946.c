/* Project Euler 946: Continued Fraction Fraction
 *
 * Sum of the first 10^8 coefficients of the continued fraction of
 * beta = (2*alpha + 3) / (3*alpha + 2), where alpha has CF digits
 * [2; 1,1,2, 1,1,1,2, 1,1,1,1,1,2, ...] with prime-length runs of 1s.
 *
 * Uses a finite-state transducer (FST) with binary lifting for fast
 * processing of long runs of 1s.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef int64_t i64;

/* ---- Floor division ---- */

static i64 floor_div(i64 n, i64 d) {
    if (d < 0) { n = -n; d = -d; }
    if (n >= 0) return n / d;
    return -((-n + d - 1) / d);
}

/* ---- State: (A, B, C, D) representing y = (A*x + B) / (C*x + D) ---- */

typedef struct { i64 A, B, C, D; } State;

static State normalize(State s) {
    if (s.C < 0 && (s.C + s.D) < 0) {
        s.A = -s.A; s.B = -s.B; s.C = -s.C; s.D = -s.D;
    }
    return s;
}

/* Consume one CF digit 'inp' of alpha, emit CF digits of beta */
static State step_transition(State s, i64 inp, i64 *outs, int *n_outs) {
    *n_outs = 0;
    s = normalize(s);

    /* Consume input: multiply by [[inp, 1], [1, 0]] */
    i64 newA = s.A * inp + s.B;
    i64 newB = s.A;
    i64 newC = s.C * inp + s.D;
    i64 newD = s.C;
    s.A = newA; s.B = newB; s.C = newC; s.D = newD;

    while (1) {
        s = normalize(s);
        i64 denom_inf = s.C;
        i64 denom_1 = s.C + s.D;
        if (denom_inf != 0 && denom_1 != 0 &&
            ((denom_inf > 0) == (denom_1 > 0))) {
            i64 q_inf = floor_div(s.A, denom_inf);
            i64 q_1 = floor_div(s.A + s.B, denom_1);
            if (q_inf == q_1) {
                i64 q = q_inf;
                outs[(*n_outs)++] = q;
                /* Extract: y' = 1/(y - q) => (C, D, A - q*C, B - q*D) */
                i64 oA = s.C, oB = s.D;
                i64 oC = s.A - q * s.C;
                i64 oD = s.B - q * s.D;
                s.A = oA; s.B = oB; s.C = oC; s.D = oD;
                continue;
            }
        }
        break;
    }
    return s;
}

/* ---- Hash table for states ---- */

#define STATE_HT_SIZE 4096
#define STATE_HT_MASK (STATE_HT_SIZE - 1)

static State states[STATE_HT_SIZE];
static int n_states = 0;
static int state_ids[STATE_HT_SIZE]; /* maps hash slot to state index */

static i64 state_hash(State s) {
    i64 h = s.A * 31 + s.B * 37 + s.C * 41 + s.D * 43;
    return h & STATE_HT_MASK;
}

static int state_find(State s) {
    i64 h = state_hash(s);
    for (int i = 0; i < STATE_HT_SIZE; i++) {
        int idx = (h + i) & STATE_HT_MASK;
        if (state_ids[idx] == -1) return -1;
        State *st = &states[state_ids[idx]];
        if (st->A == s.A && st->B == s.B && st->C == s.C && st->D == s.D)
            return state_ids[idx];
    }
    return -1;
}

static int state_insert(State s) {
    int id = n_states;
    states[id] = s;
    i64 h = state_hash(s);
    for (int i = 0; i < STATE_HT_SIZE; i++) {
        int idx = (h + i) & STATE_HT_MASK;
        if (state_ids[idx] == -1) {
            state_ids[idx] = id;
            n_states++;
            return id;
        }
    }
    return -1;
}

/* ---- Build FST and binary lifting tables ---- */

#define MAX_POW 25

static int next1[STATE_HT_SIZE], next2[STATE_HT_SIZE];
static i64 cnt1[STATE_HT_SIZE], sum1[STATE_HT_SIZE];
static i64 cnt2[STATE_HT_SIZE], sum2[STATE_HT_SIZE];
static int pow_next[MAX_POW][STATE_HT_SIZE];
static i64 pow_cnt[MAX_POW][STATE_HT_SIZE];
static i64 pow_sum[MAX_POW][STATE_HT_SIZE];

static void prepare_tables(void) {
    memset(state_ids, -1, sizeof(state_ids));
    n_states = 0;

    State init = {2, 3, 3, 2};
    int init_id = state_insert(init);

    /* BFS to collect reachable states */
    int *queue = (int *)malloc(STATE_HT_SIZE * sizeof(int));
    int qhead = 0, qtail = 0;
    queue[qtail++] = init_id;

    while (qhead < qtail) {
        int sid = queue[qhead++];
        State st = states[sid];
        for (i64 inp = 1; inp <= 2; inp++) {
            i64 outs[16];
            int n_outs;
            State nst = step_transition(st, inp, outs, &n_outs);
            int nid = state_find(nst);
            if (nid == -1) {
                nid = state_insert(nst);
                queue[qtail++] = nid;
            }
        }
    }
    free(queue);

    /* Base transitions */
    for (int sid = 0; sid < n_states; sid++) {
        State st = states[sid];
        i64 outs[16];
        int n_outs;

        State nst1 = step_transition(st, 1, outs, &n_outs);
        next1[sid] = state_find(nst1);
        cnt1[sid] = n_outs;
        sum1[sid] = 0;
        for (int j = 0; j < n_outs; j++) sum1[sid] += outs[j];

        State nst2 = step_transition(st, 2, outs, &n_outs);
        next2[sid] = state_find(nst2);
        cnt2[sid] = n_outs;
        sum2[sid] = 0;
        for (int j = 0; j < n_outs; j++) sum2[sid] += outs[j];
    }

    /* Binary lifting for repeated 1s */
    for (int sid = 0; sid < n_states; sid++) {
        pow_next[0][sid] = next1[sid];
        pow_cnt[0][sid] = cnt1[sid];
        pow_sum[0][sid] = sum1[sid];
    }
    for (int k = 1; k < MAX_POW; k++) {
        for (int sid = 0; sid < n_states; sid++) {
            int mid = pow_next[k-1][sid];
            pow_next[k][sid] = pow_next[k-1][mid];
            pow_cnt[k][sid] = pow_cnt[k-1][sid] + pow_cnt[k-1][mid];
            pow_sum[k][sid] = pow_sum[k-1][sid] + pow_sum[k-1][mid];
        }
    }
}

/* ---- Prime generator ---- */

static i64 *primes_arr = NULL;
static i64 primes_cap = 0;
static i64 n_primes = 0;
static i64 next_cand = 2;

static i64 next_prime(void) {
    if (n_primes == 0) {
        primes_arr = (i64 *)malloc(1000 * sizeof(i64));
        primes_cap = 1000;
        primes_arr[0] = 2;
        n_primes = 1;
        next_cand = 3;
        return 2;
    }
    i64 n = next_cand;
    while (1) {
        i64 r = (i64)sqrt((double)n);
        int is_p = 1;
        for (i64 i = 0; i < n_primes; i++) {
            i64 p = primes_arr[i];
            if (p > r) break;
            if (n % p == 0) { is_p = 0; break; }
        }
        if (is_p) {
            if (n_primes >= primes_cap) {
                primes_cap *= 2;
                primes_arr = (i64 *)realloc(primes_arr, primes_cap * sizeof(i64));
            }
            primes_arr[n_primes++] = n;
            next_cand = n + 2;
            return n;
        }
        n += 2;
    }
}

/* ---- Alpha digit cursor ---- */

typedef struct {
    int stage;      /* 0=initial 2, 1=ones, 2=separator 2 */
    i64 ones_left;
} AlphaCursor;

static void alpha_init(AlphaCursor *c) {
    c->stage = 0;
    c->ones_left = 0;
}

static i64 alpha_next(AlphaCursor *c) {
    if (c->stage == 0) {
        c->stage = 1;
        c->ones_left = next_prime();
        return 2;
    }
    if (c->stage == 1) {
        c->ones_left--;
        if (c->ones_left == 0) c->stage = 2;
        return 1;
    }
    /* stage 2 */
    c->stage = 1;
    c->ones_left = next_prime();
    return 2;
}

static int alpha_can_skip_ones(AlphaCursor *c) {
    return c->stage == 1 && c->ones_left > 0;
}

static void alpha_skip_ones(AlphaCursor *c, i64 k) {
    c->ones_left -= k;
    if (c->ones_left == 0) c->stage = 2;
}

/* ---- Consume ones with limit using binary lifting ---- */

static void consume_ones_with_limit(int state_id, i64 max_ones, i64 out_limit,
        int *new_state, i64 *consumed, i64 *emitted_c, i64 *emitted_s) {
    *consumed = 0; *emitted_c = 0; *emitted_s = 0;
    *new_state = state_id;

    if (max_ones <= 0 || out_limit <= 0) return;

    int sid = state_id;
    i64 cons = 0, ec = 0, es = 0;

    int bit = 63 - __builtin_clzll(max_ones);
    while (bit >= 0) {
        i64 step = (i64)1 << bit;
        if (cons + step <= max_ones) {
            i64 c = pow_cnt[bit][sid];
            if (c <= out_limit) {
                out_limit -= c;
                cons += step;
                ec += c;
                es += pow_sum[bit][sid];
                sid = pow_next[bit][sid];
            }
        }
        bit--;
    }

    *new_state = sid;
    *consumed = cons;
    *emitted_c = ec;
    *emitted_s = es;
}

/* ---- Main solver ---- */

long long p946_native(void) {
    prepare_tables();

    i64 n = 1;
    for (int i = 0; i < 8; i++) n *= 10;

    AlphaCursor cursor;
    alpha_init(&cursor);

    int state_id = 0; /* init_id is always 0 (first inserted) */
    i64 total_sum = 0;
    i64 total_cnt = 0;
    i64 cutoff = 2000;
    i64 safe_target = n - cutoff;
    if (safe_target < 0) safe_target = 0;

    /* Fast path: skip long runs of 1s */
    while (total_cnt < safe_target) {
        if (alpha_can_skip_ones(&cursor)) {
            i64 remaining_allowed = safe_target - total_cnt;
            i64 max_ones = cursor.ones_left;
            int new_state;
            i64 used, c, s;
            consume_ones_with_limit(state_id, max_ones, remaining_allowed,
                                    &new_state, &used, &c, &s);
            if (used == 0) break;
            alpha_skip_ones(&cursor, used);
            total_cnt += c;
            total_sum += s;
            state_id = new_state;
            if (used < max_ones) break;
        } else {
            i64 c = cnt2[state_id];
            if (total_cnt + c > safe_target) break;
            total_sum += sum2[state_id];
            total_cnt += c;
            state_id = next2[state_id];
            alpha_next(&cursor); /* consume the 2 */
        }
    }

    /* Finish digit by digit */
    while (total_cnt < n) {
        i64 a = alpha_next(&cursor);
        i64 outs[16];
        int n_outs;
        State st = states[state_id];
        st = step_transition(st, a, outs, &n_outs);
        state_id = state_find(st);
        for (int j = 0; j < n_outs; j++) {
            if (total_cnt >= n) break;
            total_sum += outs[j];
            total_cnt++;
        }
    }

    free(primes_arr);
    return total_sum;
}
