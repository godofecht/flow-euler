#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef long long i64;

/* floor_sum: sum_{i=0..n-1} floor((a*i + b) / m), n >= 0, m > 0 */
static i64 floor_sum(i64 n, i64 m, i64 a, i64 b) {
    i64 res = 0;
    while (1) {
        if (a >= m) {
            res += (n - 1) * n * (a / m) / 2;
            a %= m;
        }
        if (b >= m) {
            res += n * (b / m);
            b %= m;
        }
        i64 y_max = a * n + b;
        if (y_max < m) break;
        b = y_max % m;
        n = y_max / m;
        i64 tmp = m; m = a; a = tmp;
    }
    return res;
}

/* integer cube root: largest x with x^3 <= n */
static i64 icbrt(i64 n) {
    if (n <= 0) return 0;
    int bits = 0;
    i64 t = n;
    while (t) { bits++; t >>= 1; }
    i64 x = 1LL << ((bits + 2) / 3);
    while (1) {
        i64 y = (2 * x + n / (x * x)) / 3;
        if (y >= x) break;
        x = y;
    }
    while ((x + 1) * (x + 1) * (x + 1) <= n) x++;
    while (x * x * x > n) x--;
    return x;
}

/* ---- hash table for memoization (open addressing, linear probing) ---- */
#define HT_SIZE (1 << 21)
#define HT_MASK (HT_SIZE - 1)

typedef struct { i64 key; i64 val; } ht_entry;

static int ht_get(ht_entry *ht, i64 key, i64 *out) {
    i64 h = (i64)((uint64_t)key * 2654435761ULL) & (i64)HT_MASK;
    while (ht[h].key != 0) {
        if (ht[h].key == key) { *out = ht[h].val; return 1; }
        h = (h + 1) & (i64)HT_MASK;
    }
    return 0;
}

static void ht_put(ht_entry *ht, i64 key, i64 val) {
    i64 h = (i64)((uint64_t)key * 2654435761ULL) & (i64)HT_MASK;
    while (ht[h].key != 0) {
        if (ht[h].key == key) { ht[h].val = val; return; }
        h = (h + 1) & (i64)HT_MASK;
    }
    ht[h].key = key;
    ht[h].val = val;
}

/* ---- Mertens function with memoized recursion ---- */
static i64 *g_pref;
static i64 g_limit;
static ht_entry *g_ht_m;
static ht_entry *g_ht_f;

static i64 M_func(i64 n) {
    if (n <= 0) return 0;
    if (n <= g_limit) return g_pref[n];
    i64 v;
    if (ht_get(g_ht_m, n, &v)) return v;
    i64 res = 1;
    i64 l = 2;
    while (l <= n) {
        i64 q = n / l;
        i64 r = n / q;
        res -= (r - l + 1) * M_func(q);
        l = r + 1;
    }
    ht_put(g_ht_m, n, res);
    return res;
}

/* F(n) = sum_{k<=n, 3 !| k} mu(k) = M(n) + F(n/3) */
static i64 F_func(i64 n) {
    if (n <= 0) return 0;
    i64 v;
    if (ht_get(g_ht_f, n, &v)) return v;
    i64 res = M_func(n) + F_func(n / 3);
    ht_put(g_ht_f, n, res);
    return res;
}

/* Count (x,y), x>=1, y>=1, 18x+10y <= M, 3 !| y (no gcd constraint) */
static i64 count_points_nonprimitive(i64 M) {
    if (M < 28) return 0;
    i64 n = (M - 18) / 10;
    i64 b = M - 10 * n;
    i64 total = floor_sum(n, 18, 10, b);

    i64 n3 = n / 3;
    i64 b3 = M - 30 * n3;
    i64 total3 = floor_sum(n3, 18, 30, b3);

    return total - total3;
}

/* Count (x,y), x>=1, y>=1, 18x+10y <= M, gcd(x,y)=1, 3 !| y via Mobius inversion */
static i64 count_points_primitive(i64 M) {
    i64 max_d = M / 28;
    if (max_d <= 0) return 0;

    i64 limit = icbrt(max_d * max_d) + 64;

    /* linear sieve for mu */
    i64 *mu = calloc((size_t)(limit + 1), sizeof(i64));
    char *is_comp = calloc((size_t)(limit + 1), 1);
    i64 *primes = malloc((size_t)(limit + 1) * sizeof(i64));
    i64 pc = 0;
    mu[1] = 1;
    for (i64 i = 2; i <= limit; i++) {
        if (!is_comp[i]) {
            primes[pc++] = i;
            mu[i] = -1;
        }
        for (i64 j = 0; j < pc; j++) {
            i64 p = primes[j];
            i64 v = i * p;
            if (v > limit) break;
            is_comp[v] = 1;
            if (i % p == 0) {
                mu[v] = 0;
                break;
            }
            mu[v] = -mu[i];
        }
    }
    i64 *pref = calloc((size_t)(limit + 1), sizeof(i64));
    i64 s = 0;
    for (i64 i = 1; i <= limit; i++) {
        s += mu[i];
        pref[i] = s;
    }

    g_pref = pref;
    g_limit = limit;
    g_ht_m = calloc(HT_SIZE, sizeof(ht_entry));
    g_ht_f = calloc(HT_SIZE, sizeof(ht_entry));

    i64 ans = 0;
    i64 l = 1;
    while (l <= max_d) {
        i64 q = M / l;
        i64 r = M / q;
        if (r > max_d) r = max_d;
        i64 coef = F_func(r) - F_func(l - 1);
        if (coef) ans += coef * count_points_nonprimitive(q);
        l = r + 1;
    }

    free(mu);
    free(is_comp);
    free(primes);
    free(pref);
    free(g_ht_m);
    free(g_ht_f);
    return ans;
}

long long p786_native(void) {
    i64 N = 1000000000LL;
    i64 M = 3 * N + 6;
    return 2 + 4 * count_points_primitive(M);
}
