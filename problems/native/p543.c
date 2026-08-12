/* Project Euler 543 - Prime-Sum Numbers
 * S(n) = sum over x=1..n of (number of ways to write x as sum of k primes, for any k).
 * Answer = sum_{r=3..44} S(F_r) where F_r is the r-th Fibonacci number.
 *
 * Key insight: for k>=3 and n>=2k, P(n,k)=1 always. So:
 *   k=1: primePi(n)
 *   k=2: n/2-1 + primePi(n-2)-1  (Goldbach-based)
 *   k>=3: (n/2-2)*(n+1) - (n/2+1)*(n/2) + 6
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef int64_t i64;

static i64 isqrt64(i64 n) {
    if (n < 0) return 0;
    i64 r = (i64)sqrt((double)n);
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r;
}

static i64 icbrt64(i64 n) {
    if (n < 0) return 0;
    i64 r = (i64)cbrt((double)n);
    while (r * r * r > n) r--;
    while ((r + 1) * (r + 1) * (r + 1) <= n) r++;
    return r;
}

/* Sieve */
static char *sieve(i64 n) {
    char *is_comp = calloc(n + 1, 1);
    for (i64 i = 2; i * i <= n; i++) {
        if (!is_comp[i]) {
            for (i64 j = i * i; j <= n; j += i) is_comp[j] = 1;
        }
    }
    return is_comp;
}

static i64 *list_primes(i64 n, i64 *count) {
    char *ic = sieve(n);
    i64 cnt = 0;
    for (i64 i = 2; i <= n; i++) if (!ic[i]) cnt++;
    i64 *primes = malloc(cnt * sizeof(i64));
    i64 idx = 0;
    for (i64 i = 2; i <= n; i++) if (!ic[i]) primes[idx++] = i;
    free(ic);
    *count = cnt;
    return primes;
}

/* Lehmer prime counting with memoization via hash table. */
typedef struct {
    i64 key;
    i64 val;
} CacheEntry;

#define CACHE_SIZE 1048576

static CacheEntry *cache_new(void) {
    return calloc(CACHE_SIZE, sizeof(CacheEntry));
}

static i64 cache_get(CacheEntry *c, i64 key) {
    i64 h = (i64)((uint64_t)key * 2654435761ULL) & (CACHE_SIZE - 1);
    while (c[h].key != 0) {
        if (c[h].key == key) return c[h].val;
        h = (h + 1) & (CACHE_SIZE - 1);
    }
    return -1;
}

static void cache_set(CacheEntry *c, i64 key, i64 val) {
    i64 h = (i64)((uint64_t)key * 2654435761ULL) & (CACHE_SIZE - 1);
    while (c[h].key != 0) {
        if (c[h].key == key) { c[h].val = val; return; }
        h = (h + 1) & (CACHE_SIZE - 1);
    }
    c[h].key = key;
    c[h].val = val;
}

static i64 *g_primes;
static i64 *g_pi_arr;
static i64 g_limit;
static CacheEntry *g_phi_cache;
static CacheEntry *g_pi_cache;

static i64 phi(i64 x, i64 a) {
    if (a == 0) return x;
    if (a == 1) return x - x / 2;
    i64 key = x * 100 + a;
    i64 cached = cache_get(g_phi_cache, key);
    if (cached >= 0) return cached;
    i64 r = phi(x, a - 1) - phi(x / g_primes[a - 1], a - 1);
    cache_set(g_phi_cache, key, r);
    return r;
}

static i64 pi_lehmer(i64 x);

static i64 pi_lehmer(i64 x) {
    if (x <= g_limit) return g_pi_arr[x];
    i64 cached = cache_get(g_pi_cache, x);
    if (cached >= 0) return cached;

    i64 a = pi_lehmer(isqrt64(isqrt64(x)));
    i64 b = pi_lehmer(isqrt64(x));
    i64 c = pi_lehmer(icbrt64(x));

    i64 result = phi(x, a) + ((b + a - 2) * (b - a + 1)) / 2;
    for (i64 i = a + 1; i <= b; i++) {
        i64 w = x / g_primes[i - 1];
        result -= pi_lehmer(w);
        if (i <= c) {
            i64 bi = pi_lehmer(isqrt64(w));
            for (i64 j = i; j <= bi; j++) {
                result -= pi_lehmer(w / g_primes[j - 1]) - j + 1;
            }
        }
    }
    cache_set(g_pi_cache, x, result);
    return result;
}

static i64 primePi(i64 x) {
    if (x < 2) return 0;
    g_limit = isqrt64(x) + 10;
    i64 np;
    g_primes = list_primes(g_limit, &np);

    /* Build pi array up to g_limit */
    g_pi_arr = calloc(g_limit + 1, sizeof(i64));
    i64 cnt = 0;
    for (i64 i = 0; i <= g_limit; i++) {
        if (i >= 2) {
            /* Check if i is prime by scanning primes list */
        }
        g_pi_arr[i] = cnt;
        if (i >= 2) {
            /* Binary search in primes */
            i64 lo = 0, hi = np - 1, found = 0;
            while (lo <= hi) {
                i64 mid = (lo + hi) / 2;
                if (g_primes[mid] == i) { found = 1; break; }
                if (g_primes[mid] < i) lo = mid + 1;
                else hi = mid - 1;
            }
            if (found) cnt++;
            g_pi_arr[i + 1 > g_limit ? g_limit : i + 1] = cnt;
        }
    }
    /* Rebuild properly */
    cnt = 0;
    i64 pidx = 0;
    for (i64 i = 0; i <= g_limit; i++) {
        while (pidx < np && g_primes[pidx] <= i) { cnt++; pidx++; }
        g_pi_arr[i] = cnt;
    }

    g_phi_cache = cache_new();
    g_pi_cache = cache_new();

    i64 result = pi_lehmer(x);

    free(g_primes);
    free(g_pi_arr);
    free(g_phi_cache);
    free(g_pi_cache);
    return result;
}

static i64 fibonacci(i64 n) {
    i64 a = 0, b = 1;
    for (i64 i = 0; i < n; i++) {
        i64 c = a + b;
        a = b;
        b = c;
    }
    return a;
}

static i64 S(i64 n) {
    i64 res = 0;
    /* k = 1 */
    res += primePi(n);
    /* k = 2 */
    if (n >= 4) {
        res += n / 2 - 1;
        res += primePi(n - 2) - 1;
    }
    /* k >= 3 */
    if (n >= 6) {
        i64 half = n / 2;
        res += (half - 2) * (n + 1) - (half + 1) * half + 6;
    }
    return res;
}

long long p543_native(void) {
    i64 total = 0;
    for (i64 k = 3; k <= 44; k++) {
        i64 fib = fibonacci(k);
        total += S(fib);
    }
    return (long long)total;
}
