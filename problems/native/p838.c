// Project Euler 838: Not Coprime.
// f(N) = smallest positive integer not coprime to any n<=N whose LSD is 3.
// Prints ln(f(10^6)) rounded to 6 digits.
//
// Approach: force primes p%10==3 (since n=p is required) and primes p%10==7
// with p<=cbrt(N) (since p^3 ends in 3). Remaining primes p%10==7 (not forced)
// pair with primes q%10==9 via p*q<=N. Minimum weight vertex cover on this
// bipartite graph via min s-t cut (Dinic). Weights are ln(p) scaled by 1e12.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef long long i64;
typedef unsigned long long u64;

/* ---------------- Sieve ---------------- */

static char *sieve_primes(i64 n, i64 *count_out) {
    char *is_prime = (char *)malloc(n + 1);
    memset(is_prime, 1, n + 1);
    is_prime[0] = is_prime[1] = 0;
    for (i64 p = 2; p * p <= n; p++) {
        if (is_prime[p]) {
            for (i64 i = p * p; i <= n; i += p)
                is_prime[i] = 0;
        }
    }
    i64 cnt = 0;
    for (i64 i = 2; i <= n; i++) if (is_prime[i]) cnt++;
    *count_out = cnt;
    return is_prime;
}

static i64 iroot3_floor(i64 n) {
    double x = pow((double)n, 1.0 / 3.0);
    i64 r = (i64)(x + 0.5);
    while ((r + 1) * (r + 1) * (r + 1) <= n) r++;
    while (r > 0 && r * r * r > n) r--;
    return r;
}

/* ---------------- Dinic max-flow ---------------- */

typedef struct {
    i64 to;
    i64 cap;
    i64 rev; /* index of reverse edge in g[to] */
} Edge;

typedef struct {
    Edge *edges;
    i64 len;
    i64 cap;
} EdgeList;

typedef struct {
    i64 n;
    EdgeList *g;
    i64 *level;
    i64 *it;
} Dinic;

static void el_push(EdgeList *el, Edge e) {
    if (el->len == el->cap) {
        el->cap = el->cap ? el->cap * 2 : 8;
        el->edges = (Edge *)realloc(el->edges, el->cap * sizeof(Edge));
    }
    el->edges[el->len++] = e;
}

static void dinic_add_edge(Dinic *d, i64 u, i64 v, i64 cap) {
    Edge fwd = { v, cap, (i64)d->g[v].len };
    Edge rev = { u, 0,   (i64)d->g[u].len };
    el_push(&d->g[u], fwd);
    el_push(&d->g[v], rev);
}

static int dinic_bfs(Dinic *d, i64 s, i64 t) {
    for (i64 i = 0; i < d->n; i++) d->level[i] = -1;
    i64 *q = (i64 *)malloc(d->n * sizeof(i64));
    i64 head = 0, tail = 0;
    d->level[s] = 0;
    q[tail++] = s;
    while (head < tail) {
        i64 u = q[head++];
        for (i64 i = 0; i < (i64)d->g[u].len; i++) {
            Edge *e = &d->g[u].edges[i];
            if (e->cap > 0 && d->level[e->to] < 0) {
                d->level[e->to] = d->level[u] + 1;
                q[tail++] = e->to;
            }
        }
    }
    int reachable = d->level[t] >= 0;
    free(q);
    return reachable;
}

static i64 dinic_dfs(Dinic *d, i64 u, i64 t, i64 f) {
    if (u == t) return f;
    for (; d->it[u] < (i64)d->g[u].len; d->it[u]++) {
        i64 i = d->it[u];
        Edge *e = &d->g[u].edges[i];
        if (e->cap > 0 && d->level[e->to] == d->level[u] + 1) {
            i64 push = (f < e->cap) ? f : e->cap;
            i64 pushed = dinic_dfs(d, e->to, t, push);
            if (pushed > 0) {
                e->cap -= pushed;
                d->g[e->to].edges[e->rev].cap += pushed;
                return pushed;
            }
        }
    }
    return 0;
}

static i64 dinic_max_flow(Dinic *d, i64 s, i64 t) {
    i64 flow = 0;
    const i64 INF = (i64)1e18;
    while (dinic_bfs(d, s, t)) {
        for (i64 i = 0; i < d->n; i++) d->it[i] = 0;
        i64 pushed;
        while ((pushed = dinic_dfs(d, s, t, INF)) > 0)
            flow += pushed;
    }
    return flow;
}

static void dinic_reachable(Dinic *d, i64 s, char *vis) {
    memset(vis, 0, d->n);
    i64 *q = (i64 *)malloc(d->n * sizeof(i64));
    i64 head = 0, tail = 0;
    vis[s] = 1;
    q[tail++] = s;
    while (head < tail) {
        i64 u = q[head++];
        for (i64 i = 0; i < (i64)d->g[u].len; i++) {
            Edge *e = &d->g[u].edges[i];
            if (e->cap > 0 && !vis[e->to]) {
                vis[e->to] = 1;
                q[tail++] = e->to;
            }
        }
    }
    free(q);
}

static void dinic_init(Dinic *d, i64 n) {
    d->n = n;
    d->g = (EdgeList *)calloc(n, sizeof(EdgeList));
    d->level = (i64 *)malloc(n * sizeof(i64));
    d->it = (i64 *)malloc(n * sizeof(i64));
}

static void dinic_free(Dinic *d) {
    for (i64 i = 0; i < d->n; i++) free(d->g[i].edges);
    free(d->g);
    free(d->level);
    free(d->it);
}

/* ---------------- Problem logic ---------------- */

#define SCALE 1000000000000LL /* 1e12 */

double p838_native(void) {
    i64 N = 1000000;
    i64 prime_count;
    char *is_prime = sieve_primes(N, &prime_count);

    /* Collect primes into arrays by residue */
    i64 *primes = (i64 *)malloc(prime_count * sizeof(i64));
    i64 pc = 0;
    for (i64 i = 2; i <= N; i++) if (is_prime[i]) primes[pc++] = i;

    i64 cbrt = iroot3_floor(N);

    /* forced set: p%10==3, or (p%10==7 and p<=cbrt) */
    /* We'll mark forced primes in a hash-free way: they're in the primes array.
       For the bipartite graph we need left = p%10==7, not forced, p*19<=N
       and right = q%10==9, q <= N/pmin. */

    /* Build left list */
    i64 *left = (i64 *)malloc(prime_count * sizeof(i64));
    i64 nL = 0;
    for (i64 i = 0; i < pc; i++) {
        i64 p = primes[i];
        if (p % 10 == 7 && p * 19 <= N && !(p <= cbrt)) {
            /* not forced means p > cbrt (forced 7-primes are p<=cbrt) */
            left[nL++] = p;
        }
    }
    /* left is already sorted since primes are sorted */

    double ln_total = 0.0;

    /* Add forced primes' log contributions */
    for (i64 i = 0; i < pc; i++) {
        i64 p = primes[i];
        if (p % 10 == 3) {
            ln_total += log((double)p);
        } else if (p % 10 == 7 && p <= cbrt) {
            ln_total += log((double)p);
        }
    }

    if (nL == 0) {
        free(left); free(primes); free(is_prime);
        return ln_total;
    }

    i64 pmin = left[0];
    /* Build right list: q%10==9, q <= N/pmin */
    i64 *right = (i64 *)malloc(prime_count * sizeof(i64));
    i64 nR = 0;
    i64 rlim = N / pmin;
    for (i64 i = 0; i < pc; i++) {
        i64 q = primes[i];
        if (q % 10 == 9 && q <= rlim) right[nR++] = q;
    }

    if (nR == 0) {
        free(left); free(right); free(primes); free(is_prime);
        return ln_total;
    }

    /* prefix lengths: for each left[i], count right[j] <= N/left[i] */
    i64 *pref_len = (i64 *)malloc(nL * sizeof(i64));
    for (i64 i = 0; i < nL; i++) {
        i64 lim = N / left[i];
        /* binary search in right (sorted) for upper bound */
        i64 lo = 0, hi = nR;
        while (lo < hi) {
            i64 mid = (lo + hi) / 2;
            if (right[mid] <= lim) lo = mid + 1;
            else hi = mid;
        }
        pref_len[i] = lo;
    }

    /* Weights scaled by 1e12 */
    i64 *wL = (i64 *)malloc(nL * sizeof(i64));
    i64 *wR = (i64 *)malloc(nR * sizeof(i64));
    i64 total_w = 0;
    for (i64 i = 0; i < nL; i++) {
        wL[i] = (i64)(log((double)left[i]) * SCALE + 0.5);
        total_w += wL[i];
    }
    for (i64 j = 0; j < nR; j++) {
        wR[j] = (i64)(log((double)right[j]) * SCALE + 0.5);
        total_w += wR[j];
    }
    i64 INF = total_w + 1;

    /* Build flow network:
       S=0, offL=1, offR=1+nL, T=1+nL+nR */
    i64 nNodes = 2 + nL + nR;
    Dinic din;
    dinic_init(&din, nNodes);
    i64 S = 0, offL = 1, offR = 1 + nL, T = 1 + nL + nR;

    for (i64 i = 0; i < nL; i++)
        dinic_add_edge(&din, S, offL + i, wL[i]);
    for (i64 j = 0; j < nR; j++)
        dinic_add_edge(&din, offR + j, T, wR[j]);

    for (i64 i = 0; i < nL; i++) {
        i64 u = offL + i;
        i64 k = pref_len[i];
        for (i64 j = 0; j < k; j++)
            dinic_add_edge(&din, u, offR + j, INF);
    }

    dinic_max_flow(&din, S, T);

    char *vis = (char *)malloc(nNodes);
    dinic_reachable(&din, S, vis);

    /* min vertex cover = (Left not reachable) U (Right reachable) */
    for (i64 i = 0; i < nL; i++) {
        if (!vis[offL + i]) ln_total += log((double)left[i]);
    }
    for (i64 j = 0; j < nR; j++) {
        if (vis[offR + j]) ln_total += log((double)right[j]);
    }

    free(vis);
    dinic_free(&din);
    free(wL); free(wR); free(pref_len);
    free(left); free(right); free(primes); free(is_prime);
    return ln_total;
}
