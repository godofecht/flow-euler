#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef int64_t i64;
typedef uint64_t u64;
typedef uint32_t u32;

/* ----------------------------
 * Basic sieve up to sqrt(N)
 * ---------------------------- */
static int *is_prime_arr;
static u32 *primes_list;
static int plen;
static int *pi_small; /* pi_small[x] = number of primes <= x, for x <= S */
static int S; /* isqrt(N) */

static void sieve_with_pi(int limit) {
    is_prime_arr = (int *)calloc((size_t)(limit + 1), sizeof(int));
    for (int i = 0; i <= limit; i++) is_prime_arr[i] = 1;
    is_prime_arr[0] = 0;
    is_prime_arr[1] = 0;
    int r = (int)sqrt((double)limit);
    for (int i = 2; i <= r; i++) {
        if (is_prime_arr[i]) {
            for (int j = i * i; j <= limit; j += i)
                is_prime_arr[j] = 0;
        }
    }
    /* Count primes */
    int cnt = 0;
    for (int i = 0; i <= limit; i++) {
        if (is_prime_arr[i]) cnt++;
    }
    primes_list = (u32 *)calloc((size_t)cnt, sizeof(u32));
    plen = 0;
    for (int i = 2; i <= limit; i++) {
        if (is_prime_arr[i]) primes_list[plen++] = (u32)i;
    }
    /* pi_small */
    pi_small = (int *)calloc((size_t)(limit + 1), sizeof(int));
    int c = 0;
    for (int i = 0; i <= limit; i++) {
        if (is_prime_arr[i]) c++;
        pi_small[i] = c;
    }
}

/* ----------------------------
 * Prime counting table on hyperbola values
 * (optimized Lucy / hyperbola method)
 * ---------------------------- */
static i64 *g_table;
static int start_small;
static int m_table;
static i64 N_val;

static void build_prime_pi_table(i64 N, int Sval) {
    N_val = N;
    S = Sval;
    if (N / S == S) start_small = S - 1;
    else start_small = S;

    m_table = S + start_small;
    g_table = (i64 *)calloc((size_t)m_table, sizeof(i64));

    /* Large part: g[i] = floor(N/(i+1)) - 1, for i=0..S-1 */
    for (int i = 0; i < S; i++)
        g_table[i] = N / (i64)(i + 1) - 1;

    /* Small part: for v in 1..start_small, index is m-v */
    for (int v = 1; v <= start_small; v++)
        g_table[m_table - v] = (i64)(v - 1);

    /* Main Lucy updates */
    for (int pi = 0; pi < plen; pi++) {
        i64 p = primes_list[pi];
        i64 p2 = p * p;
        if (p2 > N) break;

        i64 sp = g_table[m_table - (int)(p - 1)];

        if (p2 <= (i64)S) {
            /* Update ALL large indices 0..S-1 */
            i64 j = p;
            for (int i = 0; i < S; i++) {
                if (j <= (i64)S)
                    g_table[i] -= g_table[(int)j - 1] - sp;
                else
                    g_table[i] -= g_table[m_table - (int)(N / j)] - sp;
                j += p;
            }
            /* Update small values v descending where v >= p^2 */
            if (p2 <= (i64)start_small) {
                for (int v = start_small; v >= (int)p2; v--)
                    g_table[m_table - v] -= g_table[m_table - (v / (int)p)] - sp;
            }
        } else {
            /* Only some large indices have value >= p^2 */
            i64 end = N / p2;
            if (end > (i64)S) end = (i64)S;

            i64 j = p;
            for (int i = 0; i < (int)end; i++) {
                if (j <= (i64)S)
                    g_table[i] -= g_table[(int)j - 1] - sp;
                else
                    g_table[i] -= g_table[m_table - (int)(N / j)] - sp;
                j += p;
            }
        }
    }
}

/* Prime counting function */
static i64 prime_pi(i64 x) {
    if (x <= (i64)S) return (i64)pi_small[(int)x];
    return g_table[(int)(N_val / x) - 1];
}

/* ----------------------------
 * Count squarefree numbers with exactly mleft primes
 * ---------------------------- */
/* forb: small array of forbidden primes (primes dividing the powerful part) */
static int forb_count;
static u32 forb_primes[16];

static i64 count_sqf(i64 limit, int mleft, int start_idx) {
    if (mleft == 0) return 1;
    if (start_idx >= plen) return 0;

    if (mleft == 1) {
        if (limit < (i64)primes_list[start_idx]) return 0;
        i64 base = prime_pi(limit);
        if (start_idx > 0)
            base -= (i64)pi_small[primes_list[start_idx - 1]];
        /* subtract forbidden primes in range */
        u32 start_p = primes_list[start_idx];
        for (int q = 0; q < forb_count; q++) {
            if ((i64)forb_primes[q] >= (i64)start_p && (i64)forb_primes[q] <= limit)
                base--;
        }
        return base;
    }

    if (mleft == 2) {
        i64 cnt = 0;
        for (int i = start_idx; i < plen; i++) {
            u32 p = primes_list[i];
            if ((i64)p * p > limit) break;
            /* check if p is forbidden */
            int skip = 0;
            for (int q = 0; q < forb_count; q++) {
                if (forb_primes[q] == p) { skip = 1; break; }
            }
            if (skip) continue;
            i64 lim = limit / p;
            i64 total = prime_pi(lim) - (i64)pi_small[p];
            for (int q = 0; q < forb_count; q++) {
                if ((i64)forb_primes[q] > (i64)p && (i64)forb_primes[q] <= lim)
                    total--;
            }
            cnt += total;
        }
        return cnt;
    }

    if (mleft == 3) {
        i64 cnt = 0;
        for (int i = start_idx; i < plen; i++) {
            u32 p = primes_list[i];
            if ((i64)p * p * p > limit) break;
            int skip = 0;
            for (int q = 0; q < forb_count; q++) {
                if (forb_primes[q] == p) { skip = 1; break; }
            }
            if (skip) continue;
            cnt += count_sqf(limit / p, 2, i + 1);
        }
        return cnt;
    }

    if (mleft == 4) {
        i64 cnt = 0;
        for (int i = start_idx; i < plen; i++) {
            u32 p = primes_list[i];
            i64 p4 = (i64)p * p * p * p;
            if (p4 > limit) break;
            int skip = 0;
            for (int q = 0; q < forb_count; q++) {
                if (forb_primes[q] == p) { skip = 1; break; }
            }
            if (skip) continue;
            cnt += count_sqf(limit / p, 3, i + 1);
        }
        return cnt;
    }

    return 0; /* mleft too large */
}

/* ----------------------------
 * Core solver
 * ---------------------------- */
static i64 Q[11];
static i64 N_main;

static void process_powerful(i64 a, int d) {
    i64 limit = N_main / a;
    for (int mleft = 0; mleft < 5; mleft++) {
        int tau = d << mleft;
        if (tau > 20) break;
        if (tau >= 4) {
            int k = tau / 2;
            if (k >= 2 && k <= 10) {
                Q[k] += count_sqf(limit, mleft, 0);
            }
        }
    }
}

static void dfs_powerful(int start_idx, i64 a, int d) {
    process_powerful(a, d);

    for (int i = start_idx; i < plen; i++) {
        u32 p = primes_list[i];
        if (a * (i64)p * p > N_main) break;

        i64 limit = N_main / a;
        i64 p_pow = (i64)p * p;
        int e = 2;
        while (e <= 20 && p_pow <= limit) {
            int f = (e & 1) == 0 ? e : (e + 1);
            int new_d = d * f;
            if (new_d <= 20) {
                /* Add p to forbidden primes */
                int old_fc = forb_count;
                forb_primes[forb_count++] = p;
                dfs_powerful(i + 1, a * p_pow, new_d);
                forb_count = old_fc;
            }
            e++;
            if (e > 20) break;
            p_pow *= p;
        }
    }
}

long long p861_native(void) {
    N_main = 1000000000000LL; /* 10^12 */

    int Sval = (int)sqrt((double)N_main);
    /* Adjust: isqrt */
    while ((i64)Sval * Sval > N_main) Sval--;
    while ((i64)(Sval + 1) * (Sval + 1) <= N_main) Sval++;

    sieve_with_pi(Sval);
    build_prime_pi_table(N_main, Sval);

    for (int i = 0; i < 11; i++) Q[i] = 0;
    forb_count = 0;

    dfs_powerful(0, 1, 1);

    i64 ans = 0;
    for (int k = 2; k <= 10; k++)
        ans += Q[k];

    free(is_prime_arr);
    free(primes_list);
    free(pi_small);
    free(g_table);

    return ans;
}
