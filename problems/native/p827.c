// Project Euler 827: Pythagorean Triple Occurrence
// sum_{k=1..18} Q(10^k) mod 409120391
// Uses MPFR for high-precision log2 comparison.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpfr.h>

typedef unsigned long long u64;
typedef long long i64;
typedef __int128 i128;

#define MOD 409120391LL
#define MPFR_PREC 512  // bits of precision

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

static u64 pollard_rho(u64 n) {
    if (n % 2 == 0) return 2;
    if (n % 3 == 0) return 3;
    if (n % 5 == 0) return 5;
    u64 c_vals[] = {1,3,5,7,11,13,17,19,23,29,31,37};
    for (int ci = 0; ci < 12; ci++) {
        u64 c = c_vals[ci];
        u64 x = 3, y = 3, d = 1;
        while (d == 1) {
            x = (mulmod(x, x, n) + c) % n;
            y = (mulmod(y, y, n) + c) % n;
            y = (mulmod(y, y, n) + c) % n;
            u64 diff = x > y ? x - y : y - x;
            d = gcd_u64(diff, n);
        }
        if (d != n) return d;
    }
    return n;
}

typedef struct { u64 p; int e; } Factor;
static int factorize(u64 n, Factor *out) {
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
        u64 d = pollard_rho(m);
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

static int odd_divisors(u64 n, u64 *out) {
    Factor fac[20];
    int nf = factorize(n, fac);
    Factor odd_fac[20];
    int onf = 0;
    for (int i = 0; i < nf; i++) if (fac[i].p != 2) odd_fac[onf++] = fac[i];
    int nd = 1; out[0] = 1;
    for (int i = 0; i < onf; i++) {
        int cur = nd; u64 pe = 1;
        for (int e = 1; e <= odd_fac[i].e; e++) {
            pe *= odd_fac[i].p;
            for (int j = 0; j < cur; j++) out[nd++] = out[j] * pe;
        }
    }
    for (int i = 1; i < nd; i++) { u64 key = out[i]; int j = i-1; while (j>=0 && out[j]>key) { out[j+1]=out[j]; j--; } out[j+1]=key; }
    return nd;
}

// ---- Primes by residue class ----
static int primes_p1[200];
static int primes_p3[200];
static int np1 = 0, np3 = 0;
static mpfr_t log2_p1[200], log2_p3[200];
static mpfr_t log2_2;

static void gen_primes() {
    for (int x = 2; np1 < 80 || np3 < 80; x++) {
        if (x < 2) continue;
        int is_p = 1;
        for (int d = 2; d * d <= x; d++) if (x % d == 0) { is_p = 0; break; }
        if (!is_p) continue;
        if (x % 4 == 1 && np1 < 80) {
            primes_p1[np1] = x;
            mpfr_init2(log2_p1[np1], MPFR_PREC);
            mpfr_set_ui(log2_p1[np1], x, MPFR_RNDN);
            mpfr_log2(log2_p1[np1], log2_p1[np1], MPFR_RNDN);
            np1++;
        } else if (x % 4 == 3 && np3 < 80) {
            primes_p3[np3] = x;
            mpfr_init2(log2_p3[np3], MPFR_PREC);
            mpfr_set_ui(log2_p3[np3], x, MPFR_RNDN);
            mpfr_log2(log2_p3[np3], log2_p3[np3], MPFR_RNDN);
            np3++;
        }
    }
    mpfr_init2(log2_2, MPFR_PREC);
    mpfr_set_ui(log2_2, 2, MPFR_RNDN);
    mpfr_log2(log2_2, log2_2, MPFR_RNDN);
}

// ---- Minimal representation using MPFR log2 for comparison ----
#define MAX_EXPS 50

typedef struct {
    int exps[MAX_EXPS];
    int n;
    mpfr_t log2_val;
    int initialized;
} Rep;

static void rep_init(Rep *r) { mpfr_init2(r->log2_val, MPFR_PREC); r->initialized = 1; }
static void rep_clear(Rep *r) { if (r->initialized) mpfr_clear(r->log2_val); r->initialized = 0; }
static void rep_copy(Rep *dst, const Rep *src) {
    dst->n = src->n;
    memcpy(dst->exps, src->exps, src->n * sizeof(int));
    if (!dst->initialized) mpfr_init2(dst->log2_val, MPFR_PREC);
    mpfr_set(dst->log2_val, src->log2_val, MPFR_RNDN);
    dst->initialized = 1;
}

// Memoization
typedef struct {
    u64 rem;
    int idx;
    int prev_f;
    int pset;
    Rep rep;
    int valid;
    int has_rep;
} MemoEntry;

static MemoEntry *memo_table = NULL;
static int memo_count = 0, memo_cap = 0;

static int memo_lookup(u64 rem, int idx, int prev_f, int pset, Rep *out) {
    for (int i = 0; i < memo_count; i++) {
        if (memo_table[i].rem == rem && memo_table[i].idx == idx &&
            memo_table[i].prev_f == prev_f && memo_table[i].pset == pset) {
            if (memo_table[i].valid) { rep_copy(out, &memo_table[i].rep); return 1; }
            return 0;
        }
    }
    return -1;
}

static void memo_store(u64 rem, int idx, int prev_f, int pset, Rep *rep, int valid) {
    if (memo_count >= memo_cap) {
        memo_cap = memo_cap ? memo_cap * 2 : 4096;
        memo_table = realloc(memo_table, memo_cap * sizeof(MemoEntry));
    }
    memo_table[memo_count].rem = rem;
    memo_table[memo_count].idx = idx;
    memo_table[memo_count].prev_f = prev_f;
    memo_table[memo_count].pset = pset;
    memo_table[memo_count].valid = valid;
    memo_table[memo_count].has_rep = 0;
    if (valid && rep) {
        rep_init(&memo_table[memo_count].rep);
        rep_copy(&memo_table[memo_count].rep, rep);
        memo_table[memo_count].has_rep = 1;
    }
    memo_count++;
}

static int dfs_min_rep(u64 rem, int idx, int prev_f, const int *primes,
                       const mpfr_t *logs, int nprimes, int pset, Rep *out) {
    if (rem == 1) {
        out->n = 0;
        if (!out->initialized) mpfr_init2(out->log2_val, MPFR_PREC);
        mpfr_set_ui(out->log2_val, 0, MPFR_RNDN);
        out->initialized = 1;
        return 1;
    }
    if (idx >= nprimes) return 0;

    int cached = memo_lookup(rem, idx, prev_f, pset, out);
    if (cached >= 0) return cached;

    u64 divs[20000];
    int nd = odd_divisors(rem, divs);

    Rep best;
    best.n = -1;
    best.initialized = 0;
    int found = 0;

    mpfr_t cur_log, temp;
    mpfr_init2(cur_log, MPFR_PREC);
    mpfr_init2(temp, MPFR_PREC);

    for (int i = nd - 1; i >= 0; i--) {
        u64 f = divs[i];
        if (f == 1) continue;
        if (f > (u64)prev_f) continue;
        int e = (int)((f - 1) / 2);
        if (e <= 0) continue;

        Rep sub;
        sub.initialized = 0;
        if (dfs_min_rep(rem / f, idx + 1, (int)f, primes, logs, nprimes, pset, &sub)) {
            // cur_log = logs[idx] * e + sub.log2_val
            mpfr_mul_ui(cur_log, logs[idx], e, MPFR_RNDN);
            mpfr_add(cur_log, cur_log, sub.log2_val, MPFR_RNDN);

            if (best.n < 0 || mpfr_cmp(cur_log, best.log2_val) < 0) {
                best.n = sub.n + 1;
                best.exps[0] = e;
                memcpy(best.exps + 1, sub.exps, sub.n * sizeof(int));
                if (!best.initialized) mpfr_init2(best.log2_val, MPFR_PREC);
                mpfr_set(best.log2_val, cur_log, MPFR_RNDN);
                best.initialized = 1;
            }
            found = 1;
        }
        if (sub.initialized) mpfr_clear(sub.log2_val);
    }

    mpfr_clears(cur_log, temp, NULL);

    if (found) rep_copy(out, &best);
    if (best.initialized) mpfr_clear(best.log2_val);
    memo_store(rem, idx, prev_f, pset, found ? out : NULL, found);
    return found;
}

static void min_rep_for_product(u64 P, const int *primes, const mpfr_t *logs,
                                int nprimes, int pset, Rep *out) {
    out->n = -1;
    out->initialized = 0;
    if (P == 1) {
        out->n = 0;
        mpfr_init2(out->log2_val, MPFR_PREC);
        mpfr_set_ui(out->log2_val, 0, MPFR_RNDN);
        out->initialized = 1;
        return;
    }
    dfs_min_rep(P, 0, (int)P, primes, logs, nprimes, pset, out);
}

// ---- Best rep for D ----
typedef struct {
    u64 e2;
    Rep rep3;
    mpfr_t log2_val;
    int initialized;
} DRep;

static DRep *best_d_cache = NULL;
static u64 *best_d_keys = NULL;
static int best_d_count = 0;

static void best_rep_for_D(u64 D, DRep *out) {
    if (D == 1) {
        out->e2 = 0;
        out->rep3.n = 0;
        out->rep3.initialized = 0;
        mpfr_init2(out->rep3.log2_val, MPFR_PREC);
        mpfr_set_ui(out->rep3.log2_val, 0, MPFR_RNDN);
        out->rep3.initialized = 1;
        mpfr_init2(out->log2_val, MPFR_PREC);
        mpfr_set_ui(out->log2_val, 0, MPFR_RNDN);
        out->initialized = 1;
        return;
    }

    for (int i = 0; i < best_d_count; i++) {
        if (best_d_keys[i] == D) { *out = best_d_cache[i]; return; }
    }

    u64 divs[20000];
    int nd = odd_divisors(D, divs);

    DRep best;
    best.e2 = (u64)-1;
    best.initialized = 0;
    best.rep3.initialized = 0;

    mpfr_t cur_log;
    mpfr_init2(cur_log, MPFR_PREC);

    for (int i = 0; i < nd; i++) {
        u64 a2 = divs[i];
        u64 e2 = (a2 == 1) ? 0 : (a2 + 1) / 2;
        u64 C = D / a2;
        Rep rep3;
        min_rep_for_product(C, primes_p3, log2_p3, np3, 1, &rep3);
        if (rep3.n < 0) continue;
        // cur_log = log2_2 * e2 + rep3.log2_val
        {
            mpfr_t e2_mpfr;
            mpfr_init2(e2_mpfr, MPFR_PREC);
            mpfr_set_uj(e2_mpfr, e2, MPFR_RNDN);
            mpfr_mul(cur_log, log2_2, e2_mpfr, MPFR_RNDN);
            mpfr_clear(e2_mpfr);
        }
        mpfr_add(cur_log, cur_log, rep3.log2_val, MPFR_RNDN);
        if (!best.initialized || mpfr_cmp(cur_log, best.log2_val) < 0) {
            best.e2 = e2;
            rep_copy(&best.rep3, &rep3);
            if (!best.initialized) mpfr_init2(best.log2_val, MPFR_PREC);
            mpfr_set(best.log2_val, cur_log, MPFR_RNDN);
            best.initialized = 1;
        }
        if (rep3.initialized) mpfr_clear(rep3.log2_val);
    }

    mpfr_clear(cur_log);
    *out = best;
    best_d_keys = realloc(best_d_keys, (best_d_count + 1) * sizeof(u64));
    best_d_cache = realloc(best_d_cache, (best_d_count + 1) * sizeof(DRep));
    best_d_keys[best_d_count] = D;
    best_d_cache[best_d_count] = best;
    best_d_count++;
}

// ---- Q(n) ----
typedef struct {
    Rep rep1;
    u64 e2;
    Rep rep3;
    mpfr_t log2_val;
    int initialized;
} QRep;

static QRep *q_cache = NULL;
static u64 *q_keys = NULL;
static int q_count = 0;

static void Q_rep(u64 n, QRep *out) {
    for (int i = 0; i < q_count; i++) {
        if (q_keys[i] == n) { *out = q_cache[i]; return; }
    }

    u64 S = n + 1;
    u64 divs[20000];
    int nd = odd_divisors(S, divs);

    QRep best;
    best.e2 = (u64)-1;
    best.initialized = 0;
    best.rep1.initialized = 0;
    best.rep3.initialized = 0;
    int found = 0;

    mpfr_t cur_log;
    mpfr_init2(cur_log, MPFR_PREC);

    for (int i = 0; i < nd; i++) {
        u64 B = divs[i];
        Rep rep1;
        min_rep_for_product(B, primes_p1, log2_p1, np1, 0, &rep1);
        if (rep1.n < 0) continue;
        u64 D = (2 * S) / B - 1;
        DRep drep;
        best_rep_for_D(D, &drep);
        if (!drep.initialized) continue;
        // cur_log = rep1.log2_val + drep.log2_val
        mpfr_add(cur_log, rep1.log2_val, drep.log2_val, MPFR_RNDN);
        if (!found || mpfr_cmp(cur_log, best.log2_val) < 0) {
            found = 1;
            rep_copy(&best.rep1, &rep1);
            best.e2 = drep.e2;
            rep_copy(&best.rep3, &drep.rep3);
            if (!best.initialized) mpfr_init2(best.log2_val, MPFR_PREC);
            mpfr_set(best.log2_val, cur_log, MPFR_RNDN);
            best.initialized = 1;
        }
        if (rep1.initialized) mpfr_clear(rep1.log2_val);
    }

    mpfr_clear(cur_log);
    *out = best;
    q_keys = realloc(q_keys, (q_count + 1) * sizeof(u64));
    q_cache = realloc(q_cache, (q_count + 1) * sizeof(QRep));
    q_keys[q_count] = n;
    q_cache[q_count] = best;
    q_count++;
}

static u64 Q_mod(u64 n) {
    QRep qr;
    Q_rep(n, &qr);
    u64 r = 1;
    for (int i = 0; i < qr.rep1.n; i++) {
        u64 pe = powmod(primes_p1[i], qr.rep1.exps[i], MOD);
        r = (r * pe) % MOD;
    }
    u64 p2 = powmod(2, qr.e2, MOD);
    r = (r * p2) % MOD;
    for (int i = 0; i < qr.rep3.n; i++) {
        u64 pe = powmod(primes_p3[i], qr.rep3.exps[i], MOD);
        r = (r * pe) % MOD;
    }
    return r;
}

long long p827_native(void) {
    gen_primes();
    u64 total = 0;
    for (int k = 1; k <= 18; k++) {
        u64 n = 1;
        for (int j = 0; j < k; j++) n *= 10;
        total = (total + Q_mod(n)) % MOD;
    }
    return (long long)total;
}
