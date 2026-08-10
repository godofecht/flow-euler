#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;

long long pe679_answer(void) {
    const char *alphabet = "AEFR";
    const char *patterns[4] = {"FREE", "FARE", "AREA", "REEF"};
    /* build trie */
    int next_map[64][4];
    int fail[64], out[64];
    memset(next_map, -1, sizeof(next_map));
    memset(fail, 0, sizeof(fail));
    memset(out, 0, sizeof(out));
    int nodes = 1;
    for (int pi = 0; pi < 4; pi++) {
        int v = 0;
        for (const char *p = patterns[pi]; *p; p++) {
            int c = 0; while (alphabet[c] != *p) c++;
            if (next_map[v][c] < 0) next_map[v][c] = nodes++;
            v = next_map[v][c];
        }
        out[v] |= 1 << pi;
    }
    /* fail links BFS */
    int q[64], qh = 0, qt = 0;
    for (int c = 0; c < 4; c++) {
        if (next_map[0][c] >= 0) {
            fail[next_map[0][c]] = 0;
            q[qt++] = next_map[0][c];
        } else next_map[0][c] = 0;
    }
    while (qh < qt) {
        int v = q[qh++];
        out[v] |= out[fail[v]];
        for (int c = 0; c < 4; c++) {
            int u = next_map[v][c];
            if (u >= 0) {
                fail[u] = next_map[fail[v]][c];
                q[qt++] = u;
            } else {
                next_map[v][c] = next_map[fail[v]][c];
            }
        }
    }
    int trans[64][4];
    for (int s = 0; s < nodes; s++)
        for (int c = 0; c < 4; c++)
            trans[s][c] = next_map[s][c];

    /* DP: dp[len][node][mask] — but only keep current */
    /* mask 0..15, nodes ~20 */
    int NS = nodes;
    i64 *dp = calloc(NS * 16, sizeof(i64));
    i64 *nd = calloc(NS * 16, sizeof(i64));
    dp[0 * 16 + 0] = 1;
    const int N = 30;
    for (int len = 0; len < N; len++) {
        memset(nd, 0, NS * 16 * sizeof(i64));
        for (int s = 0; s < NS; s++) {
            for (int mask = 0; mask < 16; mask++) {
                i64 cur = dp[s * 16 + mask];
                if (!cur) continue;
                for (int c = 0; c < 4; c++) {
                    int ns = trans[s][c];
                    int o = out[ns];
                    if (o & mask) continue; /* reject second match */
                    int nm = mask | o;
                    nd[ns * 16 + nm] += cur;
                }
            }
        }
        i64 *tmp = dp; dp = nd; nd = tmp;
    }
    i64 ans = 0;
    for (int s = 0; s < NS; s++) ans += dp[s * 16 + 15];
    free(dp); free(nd);
    return ans;
}
