#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;
enum { MOD = 1000000007LL };

typedef struct { int rem, lower, thr; i64 cnt; } State;

static uint64_t key_of(int rem, int lower, int thr) {
    return ((uint64_t)(uint32_t)rem << 42) | ((uint64_t)(uint32_t)lower << 21) | (uint32_t)thr;
}

/* simple open-address hash map */
enum { HCAP = 1 << 20 };
static State *ht;
static uint8_t *used;

static void hclear(void) { memset(used, 0, HCAP); }

static void hadd(int rem, int lower, int thr, i64 cnt) {
    uint64_t k = key_of(rem, lower, thr);
    uint32_t i = (uint32_t)(k * 11400714819323198485ull >> (64 - 20));
    while (used[i]) {
        if (ht[i].rem == rem && ht[i].lower == lower && ht[i].thr == thr) {
            ht[i].cnt = (ht[i].cnt + cnt) % MOD;
            return;
        }
        i = (i + 1) & (HCAP - 1);
    }
    used[i] = 1;
    ht[i] = (State){rem, lower, thr, cnt % MOD};
}

static int hdump(State *out) {
    int n = 0;
    for (int i = 0; i < HCAP; i++) if (used[i]) out[n++] = ht[i];
    return n;
}

long long pe631_answer(void) {
    const i64 n = 1000000000000000000LL;
    const int m = 40;
    ht = calloc(HCAP, sizeof(State));
    used = calloc(HCAP, 1);
    State *layer = calloc(HCAP, sizeof(State));
    State *nextb = calloc(HCAP, sizeof(State));
    hclear();
    hadd(m, 0, 0, 1);
    int ln = hdump(layer);
    i64 total = 1;
    int explicit_limit = m + 2; /* min(n, m+2) with n huge */
    for (int length = 1; length <= explicit_limit; length++) {
        hclear();
        for (int i = 0; i < ln; i++) {
            int remaining = layer[i].rem;
            int lower = layer[i].lower;
            int threshold = layer[i].thr;
            i64 count = layer[i].cnt;
            int upper = remaining + 1;
            if (upper > length) upper = length;
            for (int inv = lower; inv < upper; inv++) {
                int nrem, nlo, nthr;
                if (inv < threshold) {
                    nrem = remaining - inv;
                    nlo = inv + 1;
                    nthr = threshold + 1;
                } else {
                    nrem = remaining - inv;
                    nlo = lower;
                    nthr = inv;
                }
                hadd(nrem, nlo, nthr, count);
            }
        }
        ln = hdump(layer);
        i64 sum = 0;
        for (int i = 0; i < ln; i++) sum = (sum + layer[i].cnt) % MOD;
        total = (total + sum) % MOD;
    }
    i64 stable = 0;
    for (int i = 0; i < ln; i++) stable = (stable + layer[i].cnt) % MOD;
    i64 extra = (n - (m + 2)) % MOD;
    total = (total + extra * stable) % MOD;
    free(ht); free(used); free(layer); free(nextb);
    return total;
}
