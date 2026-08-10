#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef int64_t i64;
enum { MOD = 1000000007LL, PIVCAP = 1 << 20, W = 3 };

typedef struct { uint64_t w[W]; } Mask;

static void mask_zero(Mask *a) { a->w[0]=a->w[1]=a->w[2]=0; }
static int mask_any(const Mask *a) { return (a->w[0]|a->w[1]|a->w[2]) != 0; }
static void mask_xor(Mask *a, const Mask *b) {
    a->w[0]^=b->w[0]; a->w[1]^=b->w[1]; a->w[2]^=b->w[2];
}
static void mask_setbit(Mask *a, int bit) {
    a->w[bit >> 6] |= 1ull << (bit & 63);
}
static int mask_lsb(const Mask *a) {
    for (int i = 0; i < W; i++) if (a->w[i]) return (i << 6) + __builtin_ctzll(a->w[i]);
    return -1;
}

static i64 mod_pow(i64 a, i64 e) {
    i64 r = 1 % MOD; a %= MOD;
    while (e > 0) {
        if (e & 1) r = (i64)((__int128)r * a % MOD);
        a = (i64)((__int128)a * a % MOD);
        e >>= 1;
    }
    return r;
}

static Mask *g_basis;
static char *g_basis_used;
static int g_rank_small;

static void add_small_vector(Mask v) {
    while (mask_any(&v)) {
        int lb = mask_lsb(&v);
        if (lb < 0) return;
        if (g_basis_used[lb]) mask_xor(&v, &g_basis[lb]);
        else {
            g_basis[lb] = v;
            g_basis_used[lb] = 1;
            g_rank_small++;
            return;
        }
    }
}

typedef struct { int p; Mask mask; } Pivot;
static Pivot *g_pivots;
static uint8_t *g_pused;
static int g_pivot_count;

static int pivot_get(int p, Mask *out) {
    uint32_t i = (uint32_t)((uint64_t)(uint32_t)p * 11400714819323198485ull >> (64 - 20));
    while (g_pused[i]) {
        if (g_pivots[i].p == p) { *out = g_pivots[i].mask; return 1; }
        i = (i + 1) & (PIVCAP - 1);
    }
    return 0;
}
static void pivot_set(int p, Mask mask) {
    uint32_t i = (uint32_t)((uint64_t)(uint32_t)p * 11400714819323198485ull >> (64 - 20));
    while (g_pused[i]) {
        if (g_pivots[i].p == p) { g_pivots[i].mask = mask; return; }
        i = (i + 1) & (PIVCAP - 1);
    }
    g_pused[i] = 1;
    g_pivots[i].p = p;
    g_pivots[i].mask = mask;
    g_pivot_count++;
}

long long pe619_answer(void) {
    const int a = 1000000, b = 1234567;
    int m = b - a + 1;
    int T = (int)sqrt((double)b);
    char *is_prime = calloc(T + 1, 1);
    int *small_primes = calloc(T + 1, sizeof(int));
    int spc = 0;
    memset(is_prime, 1, T + 1);
    is_prime[0] = is_prime[1] = 0;
    for (int i = 2; (i64)i * i <= T; i++) if (is_prime[i])
        for (int j = i * i; j <= T; j += i) is_prime[j] = 0;
    for (int i = 2; i <= T; i++) if (is_prime[i]) small_primes[spc++] = i;
    int *small_bit_idx = calloc(T + 1, sizeof(int));
    for (int i = 0; i < spc; i++) small_bit_idx[small_primes[i]] = i;

    int *spf = calloc((size_t)b + 1, sizeof(int));
    int *primes = calloc(b / 5 + 10, sizeof(int));
    int pc = 0;
    for (int i = 2; i <= b; i++) {
        if (!spf[i]) { spf[i] = i; primes[pc++] = i; }
        for (int j = 0; j < pc; j++) {
            long v = (long)i * primes[j];
            if (v > b) break;
            spf[v] = primes[j];
            if (primes[j] == spf[i]) break;
        }
    }

    g_basis = calloc(spc + 1, sizeof(Mask));
    g_basis_used = calloc(spc + 1, 1);
    g_rank_small = 0;
    g_pivots = calloc(PIVCAP, sizeof(Pivot));
    g_pused = calloc(PIVCAP, 1);
    g_pivot_count = 0;

    for (int n = a; n <= b; n++) {
        int x = n;
        Mask small_mask; mask_zero(&small_mask);
        int big_prime = 0;
        while (x > 1) {
            int p = spf[x];
            int e = 0;
            while (x % p == 0) { x /= p; e++; }
            if (e & 1) {
                if (p <= T) mask_setbit(&small_mask, small_bit_idx[p]);
                else big_prime = p;
            }
        }
        if (big_prime == 0) add_small_vector(small_mask);
        else {
            Mask prev;
            if (!pivot_get(big_prime, &prev)) pivot_set(big_prime, small_mask);
            else { mask_xor(&prev, &small_mask); add_small_vector(prev); }
        }
    }

    int rank = g_pivot_count + g_rank_small;
    
    i64 ans = mod_pow(2, m - rank) - 1;
    if (ans < 0) ans += MOD;
    free(is_prime); free(small_primes); free(small_bit_idx);
    free(spf); free(primes); free(g_basis); free(g_basis_used);
    free(g_pivots); free(g_pused);
    return ans;
}
