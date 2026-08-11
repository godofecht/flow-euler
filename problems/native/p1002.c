/*
 * Project Euler 1002 - Connections II
 *
 * Given an array of 2n elements where every value appears exactly twice,
 * find the maximal number of above connections when bipartite-connecting.
 *
 * Algorithm: sweep over chords with a parity DSU and leftist heaps,
 * using a min-heap of candidate components ordered by right endpoint.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef struct {
    int right_end;
    int node;
    int root;
} Candidate;

typedef struct {
    const int *left, *right;
    int *parent, *sz, *x2p;
    int *cz, *co;
    int *hl, *hr, *hrank;
    int *ch;
} Solver;

/* ---- candidates min-heap (ordered by right_end, then node, then root) ---- */

static int cand_less(Candidate a, Candidate b) {
    if (a.right_end != b.right_end) return a.right_end < b.right_end;
    if (a.node != b.node) return a.node < b.node;
    return a.root < b.root;
}

static void cand_push(Candidate *h, int *n, Candidate c) {
    int i = (*n)++;
    h[i] = c;
    while (i > 0) {
        int p = (i - 1) / 2;
        if (cand_less(h[i], h[p])) {
            Candidate t = h[i]; h[i] = h[p]; h[p] = t;
            i = p;
        } else {
            break;
        }
    }
}

static Candidate cand_pop(Candidate *h, int *n) {
    Candidate top = h[0];
    h[0] = h[--(*n)];
    int i = 0;
    while (1) {
        int l = 2 * i + 1, r = 2 * i + 2, sm = i;
        if (l < *n && cand_less(h[l], h[sm])) sm = l;
        if (r < *n && cand_less(h[r], h[sm])) sm = r;
        if (sm == i) break;
        Candidate t = h[i]; h[i] = h[sm]; h[sm] = t;
        i = sm;
    }
    return top;
}

/* ---- parity DSU (no path compression; union by size) ---- */

static int sol_root(Solver *s, int c) {
    while (s->parent[c] != c) c = s->parent[c];
    return c;
}

static int sol_parity(Solver *s, int c, int *par_out) {
    int par = 0;
    while (s->parent[c] != c) {
        par ^= s->x2p[c];
        c = s->parent[c];
    }
    *par_out = par;
    return c;
}

/* ---- leftist heap meld (recursive; depth bounded by O(log n) rank) ---- */

static int sol_meld(Solver *s, int a, int b) {
    if (a < 0) return b;
    if (b < 0) return a;
    if (s->right[a] > s->right[b]) {
        int t = a; a = b; b = t;
    }
    s->hr[a] = sol_meld(s, s->hr[a], b);
    int lc = s->hl[a];
    int rc = s->hr[a];
    int lr = (lc < 0) ? 0 : s->hrank[lc];
    int rr = (rc < 0) ? 0 : s->hrank[rc];
    if (lr < rr) {
        s->hl[a] = rc;
        s->hr[a] = lc;
        int t = lr; lr = rr; rr = t;
    }
    s->hrank[a] = rr + 1;
    return a;
}

static int sol_pop_heap(Solver *s, int node) {
    return sol_meld(s, s->hl[node], s->hr[node]);
}

static int sol_union_opposite(Solver *s, int first, int second) {
    int fp, sp;
    int fr = sol_parity(s, first, &fp);
    int sr = sol_parity(s, second, &sp);
    if (fr == sr) {
        /* same component: parity must already be opposite */
        return fr;
    }
    if (s->sz[fr] < s->sz[sr]) {
        int t = fr; fr = sr; sr = t;
        t = fp; fp = sp; sp = t;
    }
    int shift = fp ^ sp ^ 1;
    s->parent[sr] = fr;
    s->x2p[sr] = shift;
    s->sz[fr] += s->sz[sr];
    if (shift == 0) {
        s->cz[fr] += s->cz[sr];
        s->co[fr] += s->co[sr];
    } else {
        s->cz[fr] += s->co[sr];
        s->co[fr] += s->cz[sr];
    }
    s->ch[fr] = sol_meld(s, s->ch[fr], s->ch[sr]);
    s->ch[sr] = -1;
    return fr;
}

/* ---- main sweep ---- */

static int sol_solve(Solver *s, int cc) {
    Candidate *cand = malloc((cc + 1) * sizeof(Candidate));
    int cn = 0;

    for (int chord = 0; chord < cc; chord++) {
        int le = s->left[chord];
        int re = s->right[chord];
        int current = chord;

        while (cn > 0) {
            Candidate top = cand[0];
            int rt = sol_root(s, top.root);

            if (rt != top.root) {
                cand_pop(cand, &cn);
                continue;
            }
            int chr = s->ch[rt];
            if (chr != top.node || s->right[top.node] != top.right_end) {
                cand_pop(cand, &cn);
                continue;
            }
            if (rt == sol_root(s, current)) {
                cand_pop(cand, &cn);
                continue;
            }

            if (top.right_end <= le) {
                cand_pop(cand, &cn);
                while (s->ch[rt] >= 0 && s->right[s->ch[rt]] <= le) {
                    s->ch[rt] = sol_pop_heap(s, s->ch[rt]);
                }
                int nhr = s->ch[rt];
                if (nhr >= 0) {
                    Candidate nc = { s->right[nhr], nhr, rt };
                    cand_push(cand, &cn, nc);
                }
                continue;
            }

            if (top.right_end >= re) {
                break;
            }

            cand_pop(cand, &cn);
            current = sol_union_opposite(s, chord, top.node);
        }

        current = sol_root(s, current);
        int chr = s->ch[current];
        if (chr >= 0) {
            Candidate nc = { s->right[chr], chr, current };
            cand_push(cand, &cn, nc);
        }
    }

    int answer = 0;
    for (int i = 0; i < cc; i++) {
        if (s->parent[i] == i) {
            int z = s->cz[i], o = s->co[i];
            answer += (z > o) ? z : o;
        }
    }

    free(cand);
    return answer;
}

/* ---- input parsing and entry point ---- */

long long p1002_native(void) {
    FILE *f = fopen("data/p1002_input.txt", "r");
    if (!f) return -1;

    fseek(f, 0, SEEK_END);
    long fsz = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *buf = malloc(fsz + 1);
    fread(buf, 1, fsz, f);
    buf[fsz] = '\0';
    fclose(f);

    int cap = 200000;
    int *vals = malloc(cap * sizeof(int));
    int nv = 0;
    char *p = buf;
    while (*p) {
        while (*p && (*p < '0' || *p > '9') && *p != '-') p++;
        if (!*p) break;
        int neg = 0;
        if (*p == '-') { neg = 1; p++; }
        int v = 0;
        while (*p >= '0' && *p <= '9') {
            v = v * 10 + (*p - '0');
            p++;
        }
        if (neg) v = -v;
        if (nv >= cap) {
            cap *= 2;
            vals = realloc(vals, cap * sizeof(int));
        }
        vals[nv++] = v;
    }
    free(buf);

    int maxv = 0;
    for (int i = 0; i < nv; i++)
        if (vals[i] > maxv) maxv = vals[i];

    int cc = nv / 2;

    int *left = malloc(cc * sizeof(int));
    int *right = malloc(cc * sizeof(int));
    int *oc = malloc((maxv + 1) * sizeof(int));
    memset(oc, -1, (maxv + 1) * sizeof(int));

    int nc = 0;
    for (int pos = 0; pos < nv; pos++) {
        int val = vals[pos];
        int ch = oc[val];
        if (ch < 0) {
            ch = nc;
            oc[val] = ch;
            left[ch] = pos;
            right[ch] = -1;
            nc++;
        } else {
            right[ch] = pos;
            oc[val] = -1;
        }
    }
    free(oc);
    free(vals);

    Solver s;
    s.left = left;
    s.right = right;
    s.parent = malloc(cc * sizeof(int));
    s.sz = malloc(cc * sizeof(int));
    s.x2p = malloc(cc * sizeof(int));
    s.cz = malloc(cc * sizeof(int));
    s.co = malloc(cc * sizeof(int));
    s.hl = malloc(cc * sizeof(int));
    s.hr = malloc(cc * sizeof(int));
    s.hrank = malloc(cc * sizeof(int));
    s.ch = malloc(cc * sizeof(int));

    for (int i = 0; i < cc; i++) {
        s.parent[i] = i;
        s.sz[i] = 1;
        s.x2p[i] = 0;
        s.cz[i] = 1;
        s.co[i] = 0;
        s.hl[i] = -1;
        s.hr[i] = -1;
        s.hrank[i] = 1;
        s.ch[i] = i;
    }

    int ans = sol_solve(&s, cc);

    free(s.parent); free(s.sz); free(s.x2p);
    free(s.cz); free(s.co);
    free(s.hl); free(s.hr); free(s.hrank); free(s.ch);
    free(left); free(right);

    return (long long)ans;
}
