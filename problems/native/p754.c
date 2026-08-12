#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef int32_t i32;
typedef int64_t i64;
typedef __int128 i128;

/* Problem 754: Product of Gauss Factorials.
   Compute product of GF(n) for 1 <= n <= 10^8, mod 1e9+7.
   Port of the Python reference solver. */

static const i64 MOD = 1000000007LL;
static const i64 EXP_MOD = 1000000006LL; /* MOD - 1 */
static const i64 LIMIT = 100000000LL;    /* 10^8 */

static i64 mm(i64 a, i64 b) { return (i64)((i128)a * b % MOD); }

static i64 mpow(i64 base, i64 exp) {
    i64 r = 1; base %= MOD; if (base < 0) base += MOD;
    while (exp > 0) {
        if (exp & 1) r = mm(r, base);
        base = mm(base, base);
        exp >>= 1;
    }
    return r;
}

typedef struct { i64 q, pos_prod, neg_prod, mu_sum; } Agg;

static int cmp_i64(const void *a, const void *b) {
    i64 x = *(const i64 *)a, y = *(const i64 *)b;
    return (x > y) - (x < y);
}

typedef struct { i64 key, val; } KV;

static int cmp_kv(const void *a, const void *b) {
    return (int)(*(const i64 *)a - *(const i64 *)b);
}

static KV *sf_table = NULL; /* sorted by key */
static int sf_count = 0;

static i64 sf_lookup(i64 key) {
    int lo = 0, hi = sf_count;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (sf_table[mid].key < key) lo = mid + 1;
        else hi = mid;
    }
    if (lo < sf_count && sf_table[lo].key == key) return sf_table[lo].val;
    return 1; /* shouldn't happen */
}

static void compute_superfactorials(i64 *keys, int nkeys) {
    /* sort keys, deduplicate */
    qsort(keys, nkeys, sizeof(i64), cmp_i64);
    int unique = 0;
    for (int i = 0; i < nkeys; i++) {
        if (i == 0 || keys[i] != keys[i-1]) keys[unique++] = keys[i];
    }
    nkeys = unique;

    sf_table = malloc(nkeys * sizeof(KV));
    sf_count = 0;

    int pos = 0;
    if (nkeys > 0 && keys[0] == 0) {
        sf_table[sf_count].key = 0;
        sf_table[sf_count].val = 1;
        sf_count++;
        pos = 1;
    }

    i64 stop = (nkeys > 0) ? keys[nkeys - 1] : 0;
    i64 factorial = 1, superfactorial = 1;
    for (i64 x = 1; x <= stop; x++) {
        factorial = mm(factorial, x % MOD);
        superfactorial = mm(superfactorial, factorial);
        while (pos < nkeys && keys[pos] == x) {
            sf_table[sf_count].key = x;
            sf_table[sf_count].val = superfactorial;
            sf_count++;
            pos++;
        }
    }
}

long long p754_native(void) {
    /* ---- linear sieve + aggregates ---- */
    uint8_t *mu = calloc(LIMIT + 1, 1);      /* 0=u, 1=+1, 2=-1 */
    uint8_t *comp = calloc(LIMIT + 1, 1);    /* composite flag */
    i32 *primes = malloc(6000000 * sizeof(i32));
    int pc = 0;

    if (!mu || !comp || !primes) { free(mu); free(comp); free(primes); return 0; }

    /* aggregates */
    int agg_cap = 1 << 16;
    Agg *aggs = malloc(agg_cap * sizeof(Agg));
    int agg_n = 0;

    /* first block: x=1, q=LIMIT */
    aggs[agg_n].q = LIMIT;
    aggs[agg_n].pos_prod = 1;
    aggs[agg_n].neg_prod = 1;
    aggs[agg_n].mu_sum = 1;
    agg_n++;

    mu[1] = 1; /* mu(1) = +1 */

    i64 lo = 2;
    if (lo > LIMIT) goto sieve_done;
    i64 q = LIMIT / lo;
    i64 hi = LIMIT / q;
    i64 pos_prod = 1, neg_prod = 1;
    i64 mu_sum = 0;

    for (i64 x = 2; x <= LIMIT; x++) {
        i32 mux;
        if (!comp[x]) {
            primes[pc++] = (i32)x;
            mux = 2; /* prime: mu = -1 */
            mu[x] = 2;
        } else {
            mux = mu[x];
        }

        if (mux == 1) {
            pos_prod = mm(pos_prod, x % MOD);
            mu_sum += 1;
        } else if (mux == 2) {
            neg_prod = mm(neg_prod, x % MOD);
            mu_sum -= 1;
        }

        for (int i = 0; i < pc; i++) {
            i32 p = primes[i];
            i64 y = x * p;
            if (y > LIMIT) break;
            comp[y] = 1;
            if (x % p == 0) {
                /* mu(y) = 0, already 0 from calloc */
                break;
            }
            if (mux == 1) mu[y] = 2;
            else if (mux == 2) mu[y] = 1;
            /* if mux == 0, mu[y] stays 0 */
        }

        if (x == hi) {
            if (agg_n >= agg_cap) {
                agg_cap *= 2;
                aggs = realloc(aggs, agg_cap * sizeof(Agg));
            }
            aggs[agg_n].q = q;
            aggs[agg_n].pos_prod = pos_prod;
            aggs[agg_n].neg_prod = neg_prod;
            aggs[agg_n].mu_sum = mu_sum;
            agg_n++;

            lo = hi + 1;
            if (lo > LIMIT) break;
            q = LIMIT / lo;
            hi = LIMIT / q;
            pos_prod = 1; neg_prod = 1; mu_sum = 0;
        }
    }

sieve_done:
    free(mu);
    free(comp);
    free(primes);

    /* ---- superfactorials ---- */
    i64 *keys = malloc(agg_n * sizeof(i64));
    for (int i = 0; i < agg_n; i++) keys[i] = aggs[i].q - 1;
    compute_superfactorials(keys, agg_n);
    free(keys);

    /* ---- combine ---- */
    i64 result = 1;
    for (int i = 0; i < agg_n; i++) {
        i64 aq = aggs[i].q;
        i64 exponent = (aq * (aq - 1) / 2) % EXP_MOD;
        if (exponent) {
            result = mm(result, mpow(aggs[i].pos_prod, exponent));
            i64 neg_exp = (EXP_MOD - exponent) % EXP_MOD;
            if (neg_exp)
                result = mm(result, mpow(aggs[i].neg_prod, neg_exp));
        }
        i64 sf_power = ((aggs[i].mu_sum % EXP_MOD) + EXP_MOD) % EXP_MOD;
        if (sf_power) {
            i64 sf = sf_lookup(aggs[i].q - 1);
            result = mm(result, mpow(sf, sf_power));
        }
    }

    free(aggs);
    free(sf_table);
    return (long long)result;
}
