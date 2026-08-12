// Project Euler 979: Hyperbolic frog on {7,3} heptagon tiling.
// Count length-20 closed walks on the adjacency graph.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Substitution: I -> I I II, II -> I II
// Type 0 = I, Type 1 = II

#define MAX_LAYERS 11
#define MAX_NODES 120000
#define MAX_DEG 8

static int types[MAX_LAYERS][80000];
static int type_len[MAX_LAYERS];
static int parent1[MAX_LAYERS][80000];
static int parent2[MAX_LAYERS][80000];
static int offsets[MAX_LAYERS];
static int total_nodes;

static int adj[MAX_NODES][MAX_DEG];
static int deg[MAX_NODES];

static void build_layers(int max_layer) {
    // Layer 0
    types[0][0] = 0;
    type_len[0] = 1;

    if (max_layer == 0) return;

    // Layer 1: 7 type-I tiles
    for (int i = 0; i < 7; i++) {
        types[1][i] = 0;
        parent1[1][i] = 0;
        parent2[1][i] = -1;
    }
    type_len[1] = 7;

    for (int k = 2; k <= max_layer; k++) {
        int *prev = types[k - 1];
        int m = type_len[k - 1];
        int cur_idx = 0;

        for (int j = 0; j < m; j++) {
            int t = prev[j];
            if (t == 0) {
                // I -> I I II (3 elements)
                // pos 0: I
                types[k][cur_idx] = 0;
                parent1[k][cur_idx] = j;
                parent2[k][cur_idx] = -1;
                cur_idx++;
                // pos 1: I
                types[k][cur_idx] = 0;
                parent1[k][cur_idx] = j;
                parent2[k][cur_idx] = -1;
                cur_idx++;
                // pos 2: II (trailing, boundary vertex)
                types[k][cur_idx] = 1;
                parent1[k][cur_idx] = j;
                parent2[k][cur_idx] = (j + 1) % m;
                cur_idx++;
            } else {
                // II -> I II (2 elements)
                // pos 0: I
                types[k][cur_idx] = 0;
                parent1[k][cur_idx] = j;
                parent2[k][cur_idx] = -1;
                cur_idx++;
                // pos 1: II (trailing, boundary vertex)
                types[k][cur_idx] = 1;
                parent1[k][cur_idx] = j;
                parent2[k][cur_idx] = (j + 1) % m;
                cur_idx++;
            }
        }
        type_len[k] = cur_idx;
    }
}

static void add_edge(int u, int v) {
    adj[u][deg[u]++] = v;
    adj[v][deg[v]++] = u;
}

static void build_adjacency(int max_layer) {
    memset(deg, 0, sizeof(int) * MAX_NODES);

    total_nodes = 0;
    for (int k = 0; k <= max_layer; k++) {
        offsets[k] = total_nodes;
        total_nodes += type_len[k];
    }

    int origin = offsets[0];

    // Within-layer cycle edges for layers >= 1
    for (int k = 1; k <= max_layer; k++) {
        int off = offsets[k];
        int m = type_len[k];
        for (int i = 0; i < m; i++) {
            add_edge(off + i, off + (i + 1) % m);
        }
    }

    // Layer 1 connects to origin
    if (max_layer >= 1) {
        int off1 = offsets[1];
        for (int i = 0; i < type_len[1]; i++) {
            add_edge(origin, off1 + i);
        }
    }

    // Between-layer edges from parent mappings
    for (int k = 2; k <= max_layer; k++) {
        int off = offsets[k];
        int poff = offsets[k - 1];
        for (int i = 0; i < type_len[k]; i++) {
            add_edge(off + i, poff + parent1[k][i]);
            int p = parent2[k][i];
            if (p != -1) {
                add_edge(off + i, poff + p);
            }
        }
    }
}

long long p979_native(void) {
    int n = 20;
    if (n == 0) return 1;

    int max_layer = n / 2;
    build_layers(max_layer);
    build_adjacency(max_layer);

    int N = total_nodes;
    int origin = offsets[0];

    long long *dp = calloc(N, sizeof(long long));
    long long *ndp = calloc(N, sizeof(long long));
    dp[origin] = 1;

    for (int step = 0; step < n; step++) {
        memset(ndp, 0, N * sizeof(long long));
        for (int u = 0; u < N; u++) {
            if (dp[u]) {
                for (int di = 0; di < deg[u]; di++) {
                    ndp[adj[u][di]] += dp[u];
                }
            }
        }
        long long *tmp = dp; dp = ndp; ndp = tmp;
    }

    long long result = dp[origin];
    free(dp);
    free(ndp);
    return result;
}
