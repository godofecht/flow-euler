#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef long long i64;
typedef __int128 i128;

/* Project Euler 920 - Tau Numbers
   m(k) = smallest tau number x with tau(x)=k
   M(n) = sum of all m(k) whose values do not exceed 10^n
   Compute M(16).
*/

#define LIMIT 10000000000000000LL /* 10^16 */
#define MAX_VECS 300000
#define MAX_EXP_LEN 20
#define MAX_REQ 12
#define INF128 ((i128)1 << 120)

static int primes[200];
static int num_primes;

static void sieve(void) {
    int lim = 200;
    char isp[201];
    memset(isp, 1, sizeof(isp));
    isp[0] = isp[1] = 0;
    for (int p = 2; p * p <= lim; p++)
        if (isp[p])
            for (int m = p * p; m <= lim; m += p)
                isp[m] = 0;
    num_primes = 0;
    for (int i = 2; i <= lim; i++)
        if (isp[i]) primes[num_primes++] = i;
}

static i128 ipow128(i64 base, int exp) {
    i128 result = 1, b = base;
    i128 cap = (i128)100000000000000000LL * 100; /* 10^18 */
    while (exp > 0) {
        if (exp & 1) { result *= b; if (result > cap) return cap; }
        exp >>= 1;
        if (exp > 0) { b *= b; if (b > cap) b = cap; }
    }
    return result;
}

static i128 safe_mul(i128 a, i128 b) {
    if (a == 0 || b == 0) return 0;
    i128 r = a * b;
    if (r / a != b) return INF128;
    return r;
}

/* Exponent vectors */
static int vec_len[MAX_VECS];
static int vec_data[MAX_VECS][MAX_EXP_LEN];
static int num_vecs;
static int lb_primes[20];
static int num_lb_primes;
static int dfs_exps[MAX_EXP_LEN];

static void dfs_gen(int idx, int max_e, i128 current) {
    if (idx >= num_lb_primes) return;
    int p = lb_primes[idx];
    i128 p_pow = 1;
    for (int e = 1; e <= max_e; e++) {
        p_pow *= p;
        i128 nxt = current * p_pow;
        if (nxt > LIMIT) break;
        dfs_exps[idx] = e;
        int len = idx + 1;
        vec_len[num_vecs] = len;
        memcpy(vec_data[num_vecs], dfs_exps, len * sizeof(int));
        num_vecs++;
        dfs_gen(idx + 1, e, nxt);
    }
}

/* Factorize small number */
static int fac_primes[10], fac_exps[10], fac_count;

static void factorize_small(int n) {
    fac_count = 0;
    for (int i = 0; i < num_primes && primes[i] * primes[i] <= n; i++) {
        if (n % primes[i] == 0) {
            int e = 0;
            while (n % primes[i] == 0) { n /= primes[i]; e++; }
            fac_primes[fac_count] = primes[i];
            fac_exps[fac_count] = e;
            fac_count++;
        }
    }
    if (n > 1) {
        fac_primes[fac_count] = n;
        fac_exps[fac_count] = 1;
        fac_count++;
    }
}

/* min_tau_number globals */
static int g_exps[MAX_EXP_LEN], g_r;
static i64 g_req_primes[MAX_REQ];
static int g_req_exps[MAX_REQ], g_s;
static i64 g_fillers[MAX_EXP_LEN];
static int g_num_fillers;
static i128 g_best, g_limit;

static i128 lower_bound(int mask, i128 current_prod, int next_req_idx) {
    int rem_exps[MAX_EXP_LEN], rem_count = 0;
    for (int i = 0; i < g_r; i++)
        if (!((mask >> i) & 1))
            rem_exps[rem_count++] = g_exps[i];
    /* sort descending */
    for (int i = 0; i < rem_count; i++)
        for (int j = i+1; j < rem_count; j++)
            if (rem_exps[j] > rem_exps[i]) { int t = rem_exps[i]; rem_exps[i] = rem_exps[j]; rem_exps[j] = t; }

    i64 rem_primes[MAX_EXP_LEN];
    int rp_count = 0;
    for (int i = 0; i < g_num_fillers; i++) rem_primes[rp_count++] = g_fillers[i];
    for (int i = next_req_idx; i < g_s; i++) rem_primes[rp_count++] = g_req_primes[i];
    /* sort ascending */
    for (int i = 0; i < rp_count; i++)
        for (int j = i+1; j < rp_count; j++)
            if (rem_primes[j] < rem_primes[i]) { i64 t = rem_primes[i]; rem_primes[i] = rem_primes[j]; rem_primes[j] = t; }

    i128 lb = current_prod;
    int n = rem_count < rp_count ? rem_count : rp_count;
    for (int i = 0; i < n; i++) {
        i128 pp = ipow128(rem_primes[i], rem_exps[i]);
        lb = safe_mul(lb, pp);
        if (lb >= g_best || lb >= INF128) break;
    }
    return lb;
}

static void dfs_assign(int req_idx, int mask, i128 current_prod) {
    if (current_prod >= g_best) return;
    if (req_idx == g_s) {
        int rem_exps[MAX_EXP_LEN], rem_count = 0;
        for (int i = 0; i < g_r; i++)
            if (!((mask >> i) & 1))
                rem_exps[rem_count++] = g_exps[i];
        /* sort descending */
        for (int i = 0; i < rem_count; i++)
            for (int j = i+1; j < rem_count; j++)
                if (rem_exps[j] > rem_exps[i]) { int t = rem_exps[i]; rem_exps[i] = rem_exps[j]; rem_exps[j] = t; }
        i128 total = current_prod;
        for (int i = 0; i < g_num_fillers && i < rem_count; i++) {
            total = safe_mul(total, ipow128(g_fillers[i], rem_exps[i]));
            if (total >= g_best) return;
        }
        if (total <= g_limit && total < g_best) g_best = total;
        return;
    }
    if (lower_bound(mask, current_prod, req_idx) >= g_best) return;
    i64 p = g_req_primes[req_idx];
    int need = g_req_exps[req_idx];
    int prev_e = -1;
    for (int i = 0; i < g_r; i++) {
        if ((mask >> i) & 1) continue;
        int e = g_exps[i];
        if (e < need) continue;
        if (e == prev_e) continue;
        prev_e = e;
        i128 nxt = safe_mul(current_prod, ipow128(p, e));
        if (nxt >= g_best || nxt > g_limit) continue;
        dfs_assign(req_idx + 1, mask | (1 << i), nxt);
    }
}

static i64 min_tau_number(int *exps, int r, i64 *req_pr, int *req_ex, int s) {
    if (r == 0) return 1;
    if (s > r) return -1;
    int max_exp = exps[0];
    for (int i = 0; i < s; i++)
        if (req_ex[i] > max_exp) return -1;

    /* required set */
    int req_mark[200];
    memset(req_mark, 0, sizeof(req_mark));
    for (int i = 0; i < s; i++) {
        for (int j = 0; j < num_primes; j++)
            if (primes[j] == req_pr[i]) { req_mark[j] = 1; break; }
    }
    int fillers_needed = r - s;
    g_num_fillers = 0;
    for (int j = 0; j < num_primes && g_num_fillers < fillers_needed; j++)
        if (!req_mark[j]) g_fillers[g_num_fillers++] = primes[j];

    g_r = r;
    memcpy(g_exps, exps, r * sizeof(int));
    g_s = s;
    for (int i = 0; i < s; i++) { g_req_primes[i] = req_pr[i]; g_req_exps[i] = req_ex[i]; }
    /* sort required primes descending */
    for (int i = 0; i < s; i++)
        for (int j = i+1; j < s; j++)
            if (g_req_primes[j] > g_req_primes[i]) {
                i64 tp = g_req_primes[i]; g_req_primes[i] = g_req_primes[j]; g_req_primes[j] = tp;
                int te = g_req_exps[i]; g_req_exps[i] = g_req_exps[j]; g_req_exps[j] = te;
            }
    g_limit = (i128)LIMIT;

    /* greedy upper bound */
    int avail[MAX_EXP_LEN];
    memcpy(avail, exps, r * sizeof(int));
    for (int i = 0; i < r; i++)
        for (int j = i+1; j < r; j++)
            if (avail[j] < avail[i]) { int t = avail[i]; avail[i] = avail[j]; avail[j] = t; }

    i128 gp = 1;
    int gok = 1;
    int au[MAX_EXP_LEN] = {0};
    for (int i = 0; i < s && gok; i++) {
        i64 p = g_req_primes[i];
        int need = g_req_exps[i];
        int found = -1;
        for (int j = 0; j < r; j++)
            if (!au[j] && avail[j] >= need) { found = j; break; }
        if (found < 0) { gok = 0; break; }
        au[found] = 1;
        gp = safe_mul(gp, ipow128(p, avail[found]));
        if (gp > g_limit) { gok = 0; break; }
    }
    if (gok) {
        int rem[MAX_EXP_LEN], rc = 0;
        for (int j = 0; j < r; j++) if (!au[j]) rem[rc++] = avail[j];
        for (int i = 0; i < rc; i++)
            for (int j = i+1; j < rc; j++)
                if (rem[j] > rem[i]) { int t = rem[i]; rem[i] = rem[j]; rem[j] = t; }
        for (int i = 0; i < g_num_fillers && i < rc; i++) {
            gp = safe_mul(gp, ipow128(g_fillers[i], rem[i]));
            if (gp > g_limit) { gok = 0; break; }
        }
    }
    g_best = gok ? gp : INF128;

    dfs_assign(0, 0, 1);
    if (g_best <= g_limit && g_best < INF128) return (i64)g_best;
    return -1;
}

typedef struct { i64 k, n; } kn_pair;

static int cmp_kn(const void *a, const void *b) {
    const kn_pair *pa = a, *pb = b;
    if (pa->k < pb->k) return -1;
    if (pa->k > pb->k) return 1;
    if (pa->n < pb->n) return -1;
    if (pa->n > pb->n) return 1;
    return 0;
}

long long p920_native(void) {
    sieve();

    /* generate exponent vectors */
    num_vecs = 0;
    i128 prod = 1;
    num_lb_primes = 0;
    for (int i = 0; i < num_primes; i++) {
        if (prod * primes[i] > (i128)LIMIT) break;
        prod *= primes[i];
        lb_primes[num_lb_primes++] = primes[i];
    }
    dfs_gen(0, 60, 1);

    /* collect (k, n) pairs */
    kn_pair *pairs = malloc(sizeof(kn_pair) * (num_vecs + 1));
    int np = 0;

    /* k=1, n=1 */
    pairs[np].k = 1; pairs[np].n = 1; np++;

    for (int v = 0; v < num_vecs; v++) {
        int *exps = vec_data[v];
        int r = vec_len[v];

        i64 k = 1;
        i64 req_pr[MAX_REQ];
        int req_ex[MAX_REQ];
        int s = 0;

        for (int i = 0; i < r; i++) {
            int ai1 = exps[i] + 1;
            k *= ai1;
            factorize_small(ai1);
            for (int j = 0; j < fac_count; j++) {
                int fp = fac_primes[j], fe = fac_exps[j];
                int found = -1;
                for (int t = 0; t < s; t++)
                    if (req_pr[t] == fp) { found = t; break; }
                if (found >= 0) req_ex[found] += fe;
                else { req_pr[s] = fp; req_ex[s] = fe; s++; }
            }
        }

        i64 n = min_tau_number(exps, r, req_pr, req_ex, s);
        if (n > 0 && n <= LIMIT) {
            pairs[np].k = k; pairs[np].n = n; np++;
        }
    }

    qsort(pairs, np, sizeof(kn_pair), cmp_kn);

    i128 sum = 0;
    i64 prev_k = -1;
    for (int i = 0; i < np; i++) {
        if (pairs[i].k != prev_k) {
            sum += pairs[i].n;
            prev_k = pairs[i].k;
        }
    }

    free(pairs);
    return (i64)sum;
}
