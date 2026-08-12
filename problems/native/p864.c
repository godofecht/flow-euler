#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef int64_t i64;
typedef uint64_t u64;
typedef uint32_t u32;
typedef __int128 i128;
typedef unsigned __int128 u128;

/* ----------------------------
 * Modular arithmetic helpers
 * ---------------------------- */
static u32 mod_pow_u32(u32 a, u64 e, u32 p) {
    u64 r = 1, b = a % p;
    while (e > 0) {
        if (e & 1) r = r * b % p;
        b = b * b % p;
        e >>= 1;
    }
    return (u32)r;
}

/* ----------------------------
 * Prime sieve (odd-only), primes p == 1 (mod 4)
 * ---------------------------- */
static u32 *primes_1mod4;
static int primes_1mod4_count;

static void primes_1mod4_upto(i64 limit) {
    primes_1mod4 = NULL;
    primes_1mod4_count = 0;
    if (limit < 5) return;
    int size = (int)(limit / 2 + 1);
    char *sieve = (char *)calloc((size_t)size, 1);
    for (int i = 0; i < size; i++) sieve[i] = 1;
    sieve[0] = 0;
    int r = (int)sqrt((double)limit);
    for (int p = 3; p <= r; p += 2) {
        if (sieve[p / 2]) {
            int start = (p * p) / 2;
            for (int j = start; j < size; j += p) sieve[j] = 0;
        }
    }
    int cnt = 0;
    for (int i = 1; i < size; i++) {
        if (sieve[i] && (2 * i + 1) % 4 == 1) cnt++;
    }
    primes_1mod4 = (u32 *)calloc((size_t)(cnt + 1), sizeof(u32));
    for (int i = 1; i < size; i++) {
        if (sieve[i] && (2 * i + 1) % 4 == 1)
            primes_1mod4[primes_1mod4_count++] = (u32)(2 * i + 1);
    }
    free(sieve);
}

/* ----------------------------
 * Roots of x^2 == -1 (mod p^2)
 * ---------------------------- */
static u32 sqrt_minus_one_mod_p(u32 p) {
    u64 exp = (p - 1) / 2;
    u32 g = 2;
    while (mod_pow_u32(g, exp, p) != p - 1) g++;
    return mod_pow_u32(g, (p - 1) / 4, p);
}

static void roots_minus_one_mod_p2(u32 p, u64 *r1, u64 *r2) {
    u32 r = sqrt_minus_one_mod_p(p);
    u64 p2 = (u64)p * p;
    u64 s = ((u64)r * r + 1) / p;
    u32 inv = mod_pow_u32((2 * r) % p, p - 2, p);
    u64 t = ((p - (s % p)) * inv) % p;
    u64 R = (r + t * p) % p2;
    *r1 = R;
    *r2 = (p2 - R) % p2;
}

/* ----------------------------
 * CRT combine: double the residue list
 * ---------------------------- */
static void crt_combine(i64 *out, const i64 *residues, int len,
                        i64 m, u64 a1, u64 a2, u64 p2) {
    /* Extended GCD for inverse of m mod p2 */
    u64 inv;
    {
        i128 a = m % (i128)p2;
        i128 b = (i128)p2;
        i128 old_r = a, r = b;
        i128 old_s = 1, s = 0;
        while (r != 0) {
            i128 q = old_r / r;
            i128 tmp = old_r - q * r; old_r = r; r = tmp;
            tmp = old_s - q * s; old_s = s; s = tmp;
        }
        inv = (u64)((old_s % (i128)p2 + (i128)p2) % (i128)p2);
    }
    for (int i = 0; i < len; i++) {
        i64 r = residues[i];
        u64 r_mod = (u64)(r % (i64)p2);
        u64 diff1 = (a1 + p2 - r_mod) % p2;
        u64 t1 = (u64)((u128)diff1 * inv % p2);
        out[2 * i] = (i64)((i128)r + (i128)m * (i64)t1);
        u64 diff2 = (a2 + p2 - r_mod) % p2;
        u64 t2 = (u64)((u128)diff2 * inv % p2);
        out[2 * i + 1] = (i64)((i128)r + (i128)m * (i64)t2);
    }
}

/* ----------------------------
 * Direct Mobius sum: recursive DFS
 * ---------------------------- */
static i64 dms_n, dms_D, dms_ans;
static u32 *dms_primes;
static u64 *dms_roots;
static int dms_P;

static void dms_dfs(int start_idx, i64 d, i64 mod,
                    i64 *residues, int len, int mu_sign) {
    for (int i = start_idx; i < dms_P; i++) {
        u32 p = dms_primes[i];
        i64 nd = d * (i64)p;
        if (nd > dms_D) break;

        u64 p2 = (u64)p * p;
        u64 a1 = dms_roots[2 * i];
        u64 a2 = dms_roots[2 * i + 1];
        i64 nmod = (i64)((i128)mod * (i128)p2);
        int nlen = len * 2;

        i64 *nres = (i64 *)calloc((size_t)nlen, sizeof(i64));
        crt_combine(nres, residues, len, mod, a1, a2, p2);

        /* Count A_nd(n) */
        i64 A;
        if (nmod <= dms_n) {
            i64 q = dms_n / nmod;
            i64 rem = dms_n % nmod;
            A = q * nlen;
            for (int j = 0; j < nlen; j++)
                if (nres[j] <= rem) A++;
        } else {
            A = 0;
            for (int j = 0; j < nlen; j++)
                if (nres[j] <= dms_n) A++;
        }

        dms_ans += (i64)(-mu_sign) * A;

        dms_dfs(i + 1, nd, nmod, nres, nlen, -mu_sign);
        free(nres);
    }
}

static i64 direct_mobius_sum(i64 n, i64 D) {
    primes_1mod4_upto(D);

    dms_roots = (u64 *)calloc((size_t)(primes_1mod4_count * 2), sizeof(u64));
    for (int i = 0; i < primes_1mod4_count; i++) {
        u64 r1, r2;
        roots_minus_one_mod_p2(primes_1mod4[i], &r1, &r2);
        dms_roots[2 * i] = r1;
        dms_roots[2 * i + 1] = r2;
    }

    dms_n = n;
    dms_D = D;
    dms_primes = primes_1mod4;
    dms_P = primes_1mod4_count;
    dms_ans = n; /* d=1 term */

    i64 root_res[1] = {0};
    dms_dfs(0, 1, 1, root_res, 1, 1);

    free(dms_roots);
    free(primes_1mod4);
    return dms_ans;
}

/* ----------------------------
 * Negative Pell via continued fraction
 * ---------------------------- */
static int negative_pell_fundamental(i64 D, i64 x_limit, i64 *px, i64 *py) {
    i64 a0 = (i64)sqrt((double)D);
    if (a0 * a0 == D) return 0;

    i64 m = 0, d = 1, a = a0;
    i64 p_prev = 1, p = a0;
    i64 q_prev = 0, q = 1;
    int period = 0;

    while (1) {
        m = d * a - m;
        d = (D - m * m) / d;
        a = (a0 + m) / d;

        i64 new_p = a * p + p_prev;
        i64 new_q = a * q + q_prev;
        p_prev = p; p = new_p;
        q_prev = q; q = new_q;

        period++;
        if (p_prev > x_limit) return 0;
        if (a == 2 * a0) break;
    }

    if (period % 2 == 0) return 0;
    if (p_prev * p_prev - D * q_prev * q_prev != -1) return 0;

    *px = p_prev;
    *py = q_prev;
    return 1;
}

/* ----------------------------
 * Linear sieve for Mobius up to Kmax
 * ---------------------------- */
static int *spf_kmax;
static signed char *mu_kmax;

static void linear_sieve_kmax(i64 Kmax) {
    spf_kmax = (int *)calloc((size_t)(Kmax + 1), sizeof(int));
    mu_kmax = (signed char *)calloc((size_t)(Kmax + 1), 1);
    /* Estimate prime count */
    u32 *primes = (u32 *)calloc((size_t)(Kmax / 10 + 10000), sizeof(u32));
    int pc = 0;
    mu_kmax[1] = 1;
    for (int i = 2; i <= Kmax; i++) {
        if (spf_kmax[i] == 0) {
            spf_kmax[i] = i;
            primes[pc++] = (u32)i;
            mu_kmax[i] = -1;
        }
        for (int j = 0; j < pc; j++) {
            u32 p = primes[j];
            long ip = (long)i * (long)p;
            if (ip > Kmax) break;
            spf_kmax[ip] = (int)p;
            if (i % (int)p == 0) {
                mu_kmax[ip] = 0;
                break;
            }
            mu_kmax[ip] = -mu_kmax[i];
        }
    }
    free(primes);
}

/* ----------------------------
 * Factor distinct primes of y
 * ---------------------------- */
static u32 *primes_for_fact;
static int primes_for_fact_count;

static void build_primes_for_factoring(i64 limit) {
    primes_for_fact = NULL;
    primes_for_fact_count = 0;
    if (limit < 2) return;
    int size = (int)(limit / 2 + 1);
    char *sieve = (char *)calloc((size_t)size, 1);
    for (int i = 0; i < size; i++) sieve[i] = 1;
    sieve[0] = 0;
    int r = (int)sqrt((double)limit);
    for (int p = 3; p <= r; p += 2) {
        if (sieve[p / 2]) {
            int start = (p * p) / 2;
            for (int j = start; j < size; j += p) sieve[j] = 0;
        }
    }
    int cnt = 1;
    for (int i = 1; i < size; i++) if (sieve[i]) cnt++;
    primes_for_fact = (u32 *)calloc((size_t)cnt, sizeof(u32));
    primes_for_fact[primes_for_fact_count++] = 2;
    for (int i = 1; i < size; i++)
        if (sieve[i]) primes_for_fact[primes_for_fact_count++] = (u32)(2 * i + 1);
    free(sieve);
}

static int factor_distinct_primes(i64 n, u32 *out) {
    int cnt = 0;
    i64 t = n;
    for (int i = 0; i < primes_for_fact_count; i++) {
        u32 p = primes_for_fact[i];
        if ((i64)p * p > t) break;
        if (t % p == 0) {
            out[cnt++] = p;
            while (t % p == 0) t /= p;
        }
    }
    if (t > 1) out[cnt++] = (u32)t;
    return cnt;
}

static i64 mobius_tail_sum_for_y(i64 y, i64 D) {
    u32 pf[64];
    int npf = factor_distinct_primes(y, pf);
    i64 total = 0;
    int subsets = 1 << npf;
    for (int mask = 0; mask < subsets; mask++) {
        i64 prod = 1;
        int parity = 0;
        for (int b = 0; b < npf; b++) {
            if (mask & (1 << b)) {
                prod *= pf[b];
                parity ^= 1;
            }
        }
        if (prod > D) total += parity ? -1 : 1;
    }
    return total;
}

/* ----------------------------
 * Correction via Pell
 * ---------------------------- */
static i64 correction_via_pell(i64 n, i64 D) {
    i64 Kmax = (i64)(((i128)n * n + 1) / ((i128)D * D));
    if (Kmax < 2) return 0;

    linear_sieve_kmax(Kmax);

    i64 isqrt_n = (i64)sqrt((double)n);
    build_primes_for_factoring(isqrt_n + 1);

    i64 corr = 0;

    for (i64 k = 2; k <= Kmax; k++) {
        if (mu_kmax[k] == 0) continue;

        /* Filter: if an odd prime p == 3 (mod 4) divides k, skip */
        i64 t = k;
        int ok = 1;
        while (t > 1) {
            int p = spf_kmax[t];
            t /= p;
            if (p != 2 && (p & 3) == 3) { ok = 0; break; }
        }
        if (!ok) continue;

        i64 x, y;
        if (!negative_pell_fundamental(k, n, &x, &y)) continue;

        i128 A = (i128)x * x + (i128)k * y * y;
        i128 B = (i128)2 * x * y;

        while (x <= n) {
            if (y > D)
                corr += mobius_tail_sum_for_y(y, D);
            i128 new_x = A * x + (i128)k * B * y;
            i128 new_y = B * x + A * y;
            if (new_x > (i128)n) break;
            x = (i64)new_x;
            y = (i64)new_y;
        }
    }

    free(spf_kmax);
    free(mu_kmax);
    free(primes_for_fact);
    return corr;
}

/* ----------------------------
 * Final assembly
 * ---------------------------- */
long long p864_native(void) {
    i64 N = 123567101113LL;
    i64 D = 30000000LL;

    i64 direct = direct_mobius_sum(N, D);
    i64 corr = correction_via_pell(N, D);

    return direct + corr;
}
