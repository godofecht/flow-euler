// Project Euler 673: Beds and Desks
// Count permutations commuting with two involutions B and D.
// Build graph from B/D edges, find connected components,
// classify as cycle or path, compute automorphism count.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MOD 999999937LL
#define MAXN 501

static int B[MAXN];
static int D[MAXN];
static int visited[MAXN];
static int stack_arr[MAXN];
static int comp[MAXN];

static long long fact[MAXN];

static long long mod_pow(long long base, long long exp, long long m) {
    long long r = 1, b = base % m;
    while (exp > 0) {
        if (exp & 1) r = r * b % m;
        b = b * b % m;
        exp >>= 1;
    }
    return r;
}

static void parse_file(const char *path, int *inv) {
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
    char line[256];
    while (fgets(line, sizeof(line), f)) {
        char *p = line;
        int a = 0, b = 0;
        while (*p >= '0' && *p <= '9') { a = a * 10 + (*p - '0'); p++; }
        while (*p && (*p < '0' || *p > '9')) p++;
        while (*p >= '0' && *p <= '9') { b = b * 10 + (*p - '0'); p++; }
        if (a > 0 && b > 0) {
            inv[a] = b;
            inv[b] = a;
        }
    }
    fclose(f);
}

long long p673_native(void) {
    int n = 500;

    for (int i = 0; i <= n; i++) {
        B[i] = i;
        D[i] = i;
        visited[i] = 0;
    }

    parse_file("data/p673_beds.txt", B);
    parse_file("data/p673_desks.txt", D);

    fact[0] = 1;
    for (int i = 1; i <= n; i++)
        fact[i] = fact[i - 1] * i % MOD;

    // Component types: (is_cycle, k, end_color)
    static int tc_cycle[MAXN];
    static int tc_k[MAXN];
    static int tc_color[MAXN];
    static int tc_count[MAXN];
    static long long tc_aut[MAXN];
    int num_types = 0;

    for (int start = 1; start <= n; start++) {
        if (visited[start]) continue;

        int sp = 0, cp = 0;
        stack_arr[sp++] = start;
        visited[start] = 1;

        while (sp > 0) {
            int v = stack_arr[--sp];
            comp[cp++] = v;
            int u1 = B[v], u2 = D[v];
            if (!visited[u1]) { visited[u1] = 1; stack_arr[sp++] = u1; }
            if (!visited[u2]) { visited[u2] = 1; stack_arr[sp++] = u2; }
        }

        int k = cp;
        int has_loop = 0;
        for (int i = 0; i < cp; i++) {
            int v = comp[i];
            if (B[v] == v || D[v] == v) { has_loop = 1; break; }
        }

        int is_cycle = !has_loop;
        int end_color = 0;
        long long aut;

        if (is_cycle) {
            aut = k;
        } else if (k % 2 == 1) {
            aut = 1;
        } else {
            for (int i = 0; i < cp; i++) {
                int v = comp[i];
                if (B[v] == v && D[v] != v) { end_color = 1; break; }
                if (D[v] == v && B[v] != v) { end_color = 2; break; }
            }
            aut = 2;
        }

        int found = -1;
        for (int t = 0; t < num_types; t++) {
            if (tc_cycle[t] == is_cycle && tc_k[t] == k && tc_color[t] == end_color) {
                found = t;
                break;
            }
        }
        if (found >= 0) {
            tc_count[found]++;
        } else {
            tc_cycle[num_types] = is_cycle;
            tc_k[num_types] = k;
            tc_color[num_types] = end_color;
            tc_aut[num_types] = aut;
            tc_count[num_types] = 1;
            num_types++;
        }
    }

    long long answer = 1;
    for (int t = 0; t < num_types; t++) {
        answer = answer * mod_pow(tc_aut[t], tc_count[t], MOD) % MOD;
        answer = answer * fact[tc_count[t]] % MOD;
    }

    return answer;
}
