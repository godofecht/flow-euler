// Project Euler 806: Nim on Towers of Hanoi
// Compute f(n) mod 1_000_000_007 where f(n) is the sum of indices i in the
// unique shortest 3-peg Hanoi solution for which the Nim position (disk counts
// on the three pegs) is losing for the first player.
//
// Losing indices mirror-pair to sum 2^n-1, so f(n) = k*(2^n-1)/2 mod M
// where k is the number of losing positions. k is the sum over all XOR-zero
// triples (a,b,c) with a+b+c=n of the Hanoi generating-function coefficient.
// Port of the reference Python solver.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MOD 1000000007LL
#define INV2 500000004LL

static long long *fact_arr;
static long long *inv_arr;
static long long *invfact_arr;
static long long *pow2_arr;

static long long mod_pow(long long a, long long e) {
    long long r = 1;
    a %= MOD;
    while (e) {
        if (e & 1) r = r * a % MOD;
        a = a * a % MOD;
        e >>= 1;
    }
    return r;
}

static void precompute(int nmax) {
    fact_arr = (long long *)malloc((nmax + 1) * sizeof(long long));
    inv_arr = (long long *)malloc((nmax + 1) * sizeof(long long));
    invfact_arr = (long long *)malloc((nmax + 1) * sizeof(long long));
    pow2_arr = (long long *)malloc((nmax + 1) * sizeof(long long));

    fact_arr[0] = 1;
    for (int i = 1; i <= nmax; i++) fact_arr[i] = fact_arr[i - 1] * i % MOD;

    inv_arr[1] = 1;
    for (int i = 2; i <= nmax; i++)
        inv_arr[i] = MOD - (MOD / i) * inv_arr[(int)(MOD % i)] % MOD;

    invfact_arr[nmax] = mod_pow(fact_arr[nmax], MOD - 2);
    for (int i = nmax; i >= 1; i--)
        invfact_arr[i - 1] = invfact_arr[i] * i % MOD;

    pow2_arr[0] = 1;
    for (int i = 1; i <= nmax; i++) pow2_arr[i] = pow2_arr[i - 1] * 2 % MOD;
}

// Open-addressing cache for denom_coeff keyed on (a,b,c).
typedef struct { long long key; long long val; } Entry;
static Entry *cache_table;
static int cache_cap;

static long long pack_key(int a, int b, int c) {
    // a, b, c <= 100005 < 2^17, so 51 bits total.
    return ((long long)a << 34) | ((long long)b << 17) | (long long)c;
}

static void cache_init(int cap) {
    cache_cap = cap;
    cache_table = (Entry *)malloc(cap * sizeof(Entry));
    for (int i = 0; i < cap; i++) cache_table[i].key = -1;
}

static int cache_find(long long key, long long *out) {
    int mask = cache_cap - 1;
    int h = (int)(key & mask);
    while (cache_table[h].key != -1) {
        if (cache_table[h].key == key) {
            *out = cache_table[h].val;
            return 1;
        }
        h = (h + 1) & mask;
    }
    return 0;
}

static void cache_insert(long long key, long long val) {
    int mask = cache_cap - 1;
    int h = (int)(key & mask);
    while (cache_table[h].key != -1) {
        if (cache_table[h].key == key) {
            cache_table[h].val = val;
            return;
        }
        h = (h + 1) & mask;
    }
    cache_table[h].key = key;
    cache_table[h].val = val;
}

// Coefficient of x^a y^b z^c in 1 / (1 - x^2 - y^2 - z^2 - 2xyz).
static long long denom_coeff(int a, int b, int c) {
    if (a < 0 || b < 0 || c < 0) return 0;
    if (((a ^ b) & 1) || ((a ^ c) & 1)) return 0;

    long long key = pack_key(a, b, c);
    long long cached;
    if (cache_find(key, &cached)) return cached;

    int minabc = a;
    if (b < minabc) minabc = b;
    if (c < minabc) minabc = c;

    long long ans = 0;
    int i = a & 1;
    if (i <= minabc) {
        int A = (a - i) / 2;
        int B = (b - i) / 2;
        int C = (c - i) / 2;
        int m = (a + b + c - i) / 2;

        long long term = pow2_arr[i];
        term = term * fact_arr[m] % MOD;
        term = term * invfact_arr[i] % MOD;
        term = term * invfact_arr[A] % MOD;
        term = term * invfact_arr[B] % MOD;
        term = term * invfact_arr[C] % MOD;

        while (1) {
            ans += term;
            if (ans >= MOD) ans -= MOD;

            int i2 = i + 2;
            if (i2 > minabc) break;

            // term_{i+2} / term_i = 4*A*B*C / (m*(i+1)*(i+2))
            long long ratio = (4LL * A) % MOD;
            ratio = ratio * B % MOD;
            ratio = ratio * C % MOD;
            ratio = ratio * inv_arr[m] % MOD;
            ratio = ratio * inv_arr[i + 1] % MOD;
            ratio = ratio * inv_arr[i + 2] % MOD;

            term = term * ratio % MOD;

            i = i2;
            A -= 1;
            B -= 1;
            C -= 1;
            m -= 1;
        }
    }

    cache_insert(key, ans);
    return ans;
}

// Coefficient in Fy with numerator (1+y)(1+x+z-y) = 1 + x + z + xy + yz - y^2.
static long long full_coeff(int a, int b, int c) {
    long long res = denom_coeff(a, b, c);
    res += denom_coeff(a - 1, b, c);
    if (res >= MOD) res -= MOD;
    res += denom_coeff(a, b, c - 1);
    if (res >= MOD) res -= MOD;
    res += denom_coeff(a - 1, b - 1, c);
    if (res >= MOD) res -= MOD;
    res += denom_coeff(a, b - 1, c - 1);
    if (res >= MOD) res -= MOD;
    res -= denom_coeff(a, b - 2, c);
    if (res < 0) res += MOD;
    return res;
}

// Enumerate all ordered triples (a,b,c) with a+b+c=n and a xor b xor c = 0.
// For even n, each set bit at position p>=1 gives 3 choices, so 3^popcount(n).
static int (*triples)[3];
static int n_triples;

static void gen_triples(int n) {
    if (n & 1) {
        n_triples = 0;
        return;
    }

    int bits[32];
    int nbits = 0;
    int x = n, p = 0;
    while (x) {
        if (x & 1) bits[nbits++] = p;
        x >>= 1;
        p++;
    }

    int total = 1;
    for (int bi = 0; bi < nbits; bi++)
        if (bits[bi] != 0) total *= 3;

    triples = (int (*)[3])malloc(total * sizeof(int[3]));
    triples[0][0] = 0;
    triples[0][1] = 0;
    triples[0][2] = 0;
    n_triples = 1;

    for (int bi = 0; bi < nbits; bi++) {
        int pp = bits[bi];
        if (pp == 0) continue;
        int v = 1 << (pp - 1);
        int cur = n_triples;
        for (int j = 0; j < cur; j++) {
            int a = triples[j][0], b = triples[j][1], c = triples[j][2];
            triples[j][0] = a + v;
            triples[j][1] = b + v;
            triples[j][2] = c;
            triples[cur + j][0] = a + v;
            triples[cur + j][1] = b;
            triples[cur + j][2] = c + v;
            triples[2 * cur + j][0] = a;
            triples[2 * cur + j][1] = b + v;
            triples[2 * cur + j][2] = c + v;
        }
        n_triples = 3 * cur;
    }
}

long long p806_native(void) {
    int n = 100000;

    if (n & 1) return 0;

    precompute(n + 5);
    cache_init(1 << 16);

    gen_triples(n);

    long long k = 0;
    for (int t = 0; t < n_triples; t++) {
        k += full_coeff(triples[t][0], triples[t][1], triples[t][2]);
        if (k >= MOD) k -= MOD;
    }

    long long val = (pow2_arr[n] - 1 + MOD) % MOD;
    long long result = k * val % MOD * INV2 % MOD;

    free(fact_arr);
    free(inv_arr);
    free(invfact_arr);
    free(pow2_arr);
    free(cache_table);
    free(triples);

    return result;
}
