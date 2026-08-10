#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;
enum { MOD = 1001001011LL, N = 20 };

static i64 mod_pow(i64 a, i64 e) {
    i64 r = 1 % MOD; a %= MOD;
    while (e > 0) {
        if (e & 1) r = (i64)((__int128)r * a % MOD);
        a = (i64)((__int128)a * a % MOD);
        e >>= 1;
    }
    return r;
}
static int gcd_i(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }
static int v2(int x) { int c = 0; while ((x & 1) == 0) { x >>= 1; c++; } return c; }

typedef struct {
    int lens[N], mults[N], ntypes;
    int k_cycles, tmin;
    int prefix_lt[8];
    i64 count_mod;
} PartInfo;

static PartInfo infos[2000];
static int ninfos;
static int part_buf[N];

static void emit_partition(int *p, int len, i64 *fact, i64 *invfact, i64 *inv_int, int tmax) {
    int counts[N + 1]; memset(counts, 0, sizeof(counts));
    for (int i = 0; i < len; i++) counts[p[i]]++;
    PartInfo *info = &infos[ninfos++];
    info->ntypes = 0;
    info->k_cycles = len;
    info->tmin = 100;
    for (int L = 1; L <= N; L++) if (counts[L]) {
        info->lens[info->ntypes] = L;
        info->mults[info->ntypes] = counts[L];
        info->ntypes++;
        int t = v2(L);
        if (t < info->tmin) info->tmin = t;
    }
    int count_v2[8] = {0};
    for (int L = 1; L <= N; L++) if (counts[L]) count_v2[v2(L)] += counts[L];
    int s = 0;
    for (int t = 0; t <= tmax + 1; t++) info->prefix_lt[t] = 0;
    for (int t = 0; t <= tmax; t++) {
        s += count_v2[t];
        info->prefix_lt[t + 1] = s;
    }
    i64 cm = fact[N];
    for (int L = 1; L <= N; L++) if (counts[L]) {
        int mult = counts[L];
        cm = (i64)((__int128)cm * mod_pow(inv_int[L], mult) % MOD);
        cm = (i64)((__int128)cm * invfact[mult] % MOD);
    }
    info->count_mod = cm;
}

static void gen_partitions(int rem, int maxp, int depth, i64 *fact, i64 *invfact, i64 *inv_int, int tmax) {
    if (rem == 0) { emit_partition(part_buf, depth, fact, invfact, inv_int, tmax); return; }
    for (int first = (maxp < rem ? maxp : rem); first >= 1; first--) {
        part_buf[depth] = first;
        gen_partitions(rem - first, first, depth + 1, fact, invfact, inv_int, tmax);
    }
}

long long pe626_answer(void) {
    i64 fact[N + 1], invfact[N + 1], inv_int[N + 1];
    fact[0] = 1;
    for (int i = 1; i <= N; i++) fact[i] = (i64)((__int128)fact[i - 1] * i % MOD);
    invfact[N] = mod_pow(fact[N], MOD - 2);
    for (int i = N; i > 0; i--) invfact[i - 1] = (i64)((__int128)invfact[i] * i % MOD);
    for (int i = 1; i <= N; i++) inv_int[i] = mod_pow(i, MOD - 2);
    int tmax = 0; for (int p = 1; p * 2 <= N; p *= 2) tmax++;
    ninfos = 0;
    gen_partitions(N, N, 0, fact, invfact, inv_int, tmax);

    int gcd_tab[N + 1][N + 1];
    for (int i = 1; i <= N; i++)
        for (int j = 1; j <= N; j++)
            gcd_tab[i][j] = gcd_i(i, j);
    i64 pow2[N * N + 1];
    pow2[0] = 1;
    for (int i = 1; i <= N * N; i++) pow2[i] = (pow2[i - 1] * 2) % MOD;

    i64 total = 0;
    for (int i = 0; i < ninfos; i++) {
        PartInfo *pr = &infos[i];
        for (int j = 0; j < ninfos; j++) {
            PartInfo *pc = &infos[j];
            int cycles = 0;
            for (int a = 0; a < pr->ntypes; a++)
                for (int b = 0; b < pc->ntypes; b++)
                    cycles += pr->mults[a] * pc->mults[b] * gcd_tab[pr->lens[a]][pc->lens[b]];
            int d;
            if (pr->tmin < pc->tmin) d = pr->prefix_lt[pc->tmin];
            else if (pc->tmin < pr->tmin) d = pc->prefix_lt[pr->tmin];
            else d = 1;
            int e = cycles - pr->k_cycles - pc->k_cycles + d;
            i64 term = (i64)((__int128)pr->count_mod * pc->count_mod % MOD);
            term = (i64)((__int128)term * pow2[e] % MOD);
            total = (total + term) % MOD;
        }
    }
    i64 inv_fact_n = mod_pow(fact[N], MOD - 2);
    i64 inv_den = (i64)((__int128)inv_fact_n * inv_fact_n % MOD);
    return (i64)((__int128)total * inv_den % MOD);
}
