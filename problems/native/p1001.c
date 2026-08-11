// Project Euler 1001: Connections I
// Connectivity number of a chord diagram = number of independent sets
// in the chord-intersection graph, mod 1003443221.
// Port of cirosantilli's Python solver to C.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MOD 1003443221LL

static int *left_end;
static int *right_end;
static int *chord_at;

// DSU
static int *parent;
static int *sz;

// Leftist heap (one node per chord)
static int *hleft;
static int *hright;
static int *hrank;
static int *comp_heap;

static int find(int x) {
    int root = x;
    while (parent[root] != root) root = parent[root];
    while (parent[x] != root) {
        int nxt = parent[x];
        parent[x] = root;
        x = nxt;
    }
    return root;
}

static int meld(int a, int b) {
    if (a < 0) return b;
    if (b < 0) return a;
    if (right_end[a] > right_end[b]) { int t = a; a = b; b = t; }
    hright[a] = meld(hright[a], b);
    int lc = hleft[a], rc = hright[a];
    int lr = (lc < 0) ? 0 : hrank[lc];
    int rr = (rc < 0) ? 0 : hrank[rc];
    if (lr < rr) { hleft[a] = rc; hright[a] = lc; }
    hrank[a] = ((lr < rr) ? lr : rr) + 1;
    return a;
}

static int pop_heap(int node) {
    return meld(hleft[node], hright[node]);
}

static int union_dsu(int a, int b) {
    a = find(a); b = find(b);
    if (a == b) return a;
    if (sz[a] < sz[b]) { int t = a; a = b; b = t; }
    parent[b] = a;
    sz[a] += sz[b];
    comp_heap[a] = meld(comp_heap[a], comp_heap[b]);
    comp_heap[b] = -1;
    return a;
}

// Min-heap of (right, node, root) tuples
typedef struct { int right; int node; int root; } Cand;
static Cand *cand_heap;
static int cand_size;

static void cand_push(int right, int node, int root) {
    int i = cand_size++;
    cand_heap[i].right = right;
    cand_heap[i].node = node;
    cand_heap[i].root = root;
    while (i > 0) {
        int p = (i - 1) / 2;
        if (cand_heap[p].right <= cand_heap[i].right) break;
        Cand t = cand_heap[p]; cand_heap[p] = cand_heap[i]; cand_heap[i] = t;
        i = p;
    }
}

static Cand cand_pop(void) {
    Cand top = cand_heap[0];
    cand_heap[0] = cand_heap[--cand_size];
    int i = 0;
    while (1) {
        int l = 2*i+1, r = 2*i+2, s = i;
        if (l < cand_size && cand_heap[l].right < cand_heap[s].right) s = l;
        if (r < cand_size && cand_heap[r].right < cand_heap[s].right) s = r;
        if (s == i) break;
        Cand t = cand_heap[s]; cand_heap[s] = cand_heap[i]; cand_heap[i] = t;
        i = s;
    }
    return top;
}

// Parse chords from values array
static void parse_chords(int *values, int n) {
    int chord_count = n / 2;
    left_end = (int*)malloc(chord_count * sizeof(int));
    right_end = (int*)malloc(chord_count * sizeof(int));
    chord_at = (int*)malloc(n * sizeof(int));

    // Simple hash map for open chords (value -> chord index)
    // Values can be up to 20000, use direct array
    int max_val = 0;
    for (int i = 0; i < n; i++) if (values[i] > max_val) max_val = values[i];
    int *open_chord = (int*)malloc((max_val + 1) * sizeof(int));
    memset(open_chord, -1, (max_val + 1) * sizeof(int));

    int nc = 0;
    for (int pos = 0; pos < n; pos++) {
        int v = values[pos];
        int c = open_chord[v];
        if (c < 0) {
            c = nc++;
            open_chord[v] = c;
            left_end[c] = pos;
            right_end[c] = -1;
        } else {
            right_end[c] = pos;
            open_chord[v] = -1;
        }
        chord_at[pos] = c;
    }
    free(open_chord);
}

// Find crossing components
static void find_components(int chord_count) {
    parent = (int*)malloc(chord_count * sizeof(int));
    sz = (int*)malloc(chord_count * sizeof(int));
    hleft = (int*)malloc(chord_count * sizeof(int));
    hright = (int*)malloc(chord_count * sizeof(int));
    hrank = (int*)malloc(chord_count * sizeof(int));
    comp_heap = (int*)malloc(chord_count * sizeof(int));
    cand_heap = (Cand*)malloc((chord_count + 1) * sizeof(Cand));
    cand_size = 0;

    for (int i = 0; i < chord_count; i++) {
        parent[i] = i;
        sz[i] = 1;
        hleft[i] = -1;
        hright[i] = -1;
        hrank[i] = 1;
        comp_heap[i] = i;
    }

    for (int chord = 0; chord < chord_count; chord++) {
        int le = left_end[chord];
        int re = right_end[chord];
        int current = chord;

        while (cand_size > 0) {
            Cand top = cand_heap[0];
            int root = find(top.root);

            if (root != top.root) {
                cand_pop();
                continue;
            }
            int hr = comp_heap[root];
            if (hr != top.node || right_end[top.node] != top.right) {
                cand_pop();
                continue;
            }
            if (root == find(current)) {
                cand_pop();
                continue;
            }

            if (top.right <= le) {
                cand_pop();
                while (comp_heap[root] >= 0 && right_end[comp_heap[root]] <= le) {
                    comp_heap[root] = pop_heap(comp_heap[root]);
                }
                hr = comp_heap[root];
                if (hr >= 0) {
                    cand_push(right_end[hr], hr, root);
                }
                continue;
            }

            if (top.right >= re) break;

            cand_pop();
            current = union_dsu(current, root);
        }

        current = find(current);
        int hr = comp_heap[current];
        cand_push(right_end[hr], hr, current);
    }

    // Path compress all
    for (int i = 0; i < chord_count; i++) find(i);
}

// Permutation cut: find rotation where first half has each chord once
static int permutation_cut(int *word, int wc_len) {
    int chord_count = wc_len / 2;
    // Use a frequency array. Values in word are chord IDs (0..chord_count-1)
    int *freq = (int*)calloc(chord_count, sizeof(int));
    int repeated = 0;

    // doubled = word + word
    for (int i = 0; i < chord_count; i++) {
        int c = word[i];
        if (freq[c] == 1) repeated++;
        freq[c]++;
    }

    for (int start = 0; start < chord_count; start++) {
        if (repeated == 0) {
            free(freq);
            return start;
        }
        int outgoing = word[start];
        if (freq[outgoing] == 2) repeated--;
        freq[outgoing]--;
        int incoming = word[start + chord_count];
        if (freq[incoming] == 1) repeated++;
        freq[incoming]++;
    }

    free(freq);
    return -1;
}

// Count via Fenwick tree for permutation diagrams
static long long count_perm_diagram(int *word, int wc_len, int start) {
    int cc = wc_len / 2;
    long long *fenwick = (long long*)calloc(cc + 1, sizeof(long long));
    long long sub_count = 0;

    // Build second_position: chord -> position in second half
    int *second_pos = (int*)malloc(cc * sizeof(int));
    for (int i = 0; i < cc; i++) {
        int chord = word[start + cc + i];
        second_pos[chord] = i;
    }

    for (int i = 0; i < cc; i++) {
        int chord = word[start + i];
        int pos = second_pos[chord];

        long long prefix = 0;
        int idx = pos + 1;
        while (idx > 0) {
            prefix = (prefix + fenwick[idx]) % MOD;
            idx -= idx & (-idx);
        }

        long long ending = (1 + sub_count - prefix) % MOD;
        if (ending < 0) ending += MOD;
        sub_count = (sub_count + ending) % MOD;

        idx = pos + 1;
        while (idx <= cc) {
            fenwick[idx] = (fenwick[idx] + ending) % MOD;
            idx += idx & (-idx);
        }
    }

    free(fenwick);
    free(second_pos);
    return (sub_count + 1) % MOD;
}

// Count via interval DP
static long long count_interval_dp(int *word, int wc_len) {
    int n = wc_len;
    int *mate = (int*)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) mate[i] = -1;

    // Find mates using a simple approach (chord IDs are 0..cc-1)
    int cc = n / 2;
    int *first_pos = (int*)malloc(cc * sizeof(int));
    for (int i = 0; i < cc; i++) first_pos[i] = -1;

    for (int pos = 0; pos < n; pos++) {
        int c = word[pos];
        if (first_pos[c] < 0) {
            first_pos[c] = pos;
        } else {
            mate[first_pos[c]] = pos;
            mate[pos] = first_pos[c];
        }
    }

    long long *inside = (long long*)malloc(n * sizeof(long long));
    long long *current = (long long*)malloc((n + 1) * sizeof(long long));
    for (int i = 0; i < n; i++) inside[i] = 1;

    for (int boundary = 0; boundary < n; boundary++) {
        current[boundary + 1] = 1;
        for (int pos = boundary; pos >= 0; pos--) {
            int other = mate[pos];
            long long val = current[pos + 1];
            if (pos < other && other <= boundary) {
                val = (val + inside[pos] * current[other + 1]) % MOD;
            }
            current[pos] = val;
        }
        int next_pos = boundary + 1;
        if (next_pos < n && mate[next_pos] < next_pos) {
            int le = mate[next_pos];
            inside[le] = current[le + 1];
        }
    }

    long long result = current[0];
    free(mate);
    free(first_pos);
    free(inside);
    free(current);
    return result;
}

static long long count_component(int *word, int wc_len) {
    int start = permutation_cut(word, wc_len);
    if (start >= 0) return count_perm_diagram(word, wc_len, start);
    return count_interval_dp(word, wc_len);
}

long long p1001_native(void) {
    // Read input file
    FILE *f = fopen("data/p1001_input.txt", "r");
    if (!f) { fprintf(stderr, "cannot open data/p1001_input.txt\n"); return -1; }

    // Read entire file
    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *buf = (char*)malloc(fsize + 1);
    fread(buf, 1, fsize, f);
    buf[fsize] = '\0';
    fclose(f);

    // Parse comma-separated integers
    int cap = 40001;
    int *values = (int*)malloc(cap * sizeof(int));
    int n = 0;
    char *p = buf;
    while (*p) {
        // skip whitespace
        while (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t') p++;
        if (*p == '\0') break;
        int val = 0;
        while (*p >= '0' && *p <= '9') {
            val = val * 10 + (*p - '0');
            p++;
        }
        values[n++] = val;
        while (*p && *p != ',') p++;
        if (*p == ',') p++;
    }
    free(buf);

    int chord_count = n / 2;

    parse_chords(values, n);
    find_components(chord_count);

    // Group chords by root
    // Build words for each component
    int *comp_word_count = (int*)calloc(chord_count, sizeof(int));
    for (int i = 0; i < n; i++) {
        int r = find(chord_at[i]);
        comp_word_count[r]++;
    }

    // Allocate word arrays
    int **comp_words = (int**)calloc(chord_count, sizeof(int*));
    for (int i = 0; i < chord_count; i++) {
        if (comp_word_count[i] > 0) {
            comp_words[i] = (int*)malloc(comp_word_count[i] * sizeof(int));
        }
    }
    int *comp_idx = (int*)calloc(chord_count, sizeof(int));
    for (int i = 0; i < n; i++) {
        int r = find(chord_at[i]);
        comp_words[r][comp_idx[r]++] = chord_at[i];
    }

    // Multiply component counts
    long long answer = 1;
    for (int i = 0; i < chord_count; i++) {
        if (comp_word_count[i] > 0) {
            long long c = count_component(comp_words[i], comp_word_count[i]);
            answer = answer * c % MOD;
        }
    }

    // Cleanup
    for (int i = 0; i < chord_count; i++) {
        if (comp_words[i]) free(comp_words[i]);
    }
    free(comp_words);
    free(comp_idx);
    free(comp_word_count);
    free(left_end);
    free(right_end);
    free(chord_at);
    free(parent);
    free(sz);
    free(hleft);
    free(hright);
    free(hrank);
    free(comp_heap);
    free(cand_heap);
    free(values);

    return answer;
}
