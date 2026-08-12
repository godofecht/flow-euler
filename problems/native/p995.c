// Project Euler 995: S(p) product over primes < 20000.
// The final product has 536281 digits, so we use MPFR for log10 tracking
// and GMP mpz for exact product computation, then format in scientific notation.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>

#define LIMIT 20000
#define PRIME_SEARCH_LIMIT 2000000

// Sieve and prime list.
static char *is_prime_arr;
static int *primes;
static int num_primes;

static void sieve(int n) {
    is_prime_arr = calloc(n + 1, 1);
    memset(is_prime_arr, 1, n + 1);
    is_prime_arr[0] = is_prime_arr[1] = 0;
    for (int i = 2; (long long)i * i <= n; i++) {
        if (is_prime_arr[i]) {
            for (long long j = (long long)i * i; j <= n; j += i)
                is_prime_arr[j] = 0;
        }
    }
    primes = malloc((n / 8 + 1000) * sizeof(int));
    num_primes = 0;
    for (int i = 2; i <= n; i++)
        if (is_prime_arr[i]) primes[num_primes++] = i;
}

// Factor n using the prime list.
typedef struct { int p; int e; } Factor;
static int factorize(int n, Factor *out) {
    int count = 0;
    int t = n;
    for (int i = 0; i < num_primes; i++) {
        int p = primes[i];
        if ((long long)p * p > t) break;
        if (t % p == 0) {
            int e = 0;
            while (t % p == 0) { t /= p; e++; }
            out[count].p = p; out[count].e = e; count++;
        }
    }
    if (t > 1) { out[count].p = t; out[count].e = 1; count++; }
    return count;
}

// Divisors from factorization.
static int divisors_from_factors(const Factor *factors, int nf, int *divs) {
    int ndivs = 1;
    divs[0] = 1;
    for (int i = 0; i < nf; i++) {
        int p = factors[i].p, e = factors[i].e;
        int old_ndivs = ndivs;
        int power = 1;
        for (int j = 0; j < e; j++) {
            power *= p;
            for (int k = 0; k < old_ndivs; k++)
                divs[ndivs++] = divs[k] * power;
        }
    }
    // Sort divisors.
    // Simple sort (small arrays).
    for (int i = 1; i < ndivs; i++) {
        int key = divs[i];
        int j = i - 1;
        while (j >= 0 && divs[j] > key) { divs[j+1] = divs[j]; j--; }
        divs[j+1] = key;
    }
    return ndivs;
}

static int gcd_int(int a, int b) {
    while (b) { int t = a % b; a = b; b = t; }
    return a;
}

// Primitive root modulo prime p.
static int primitive_root(int p, const Factor *factors_pm1, int nf) {
    if (p == 2) return 1;
    int m = p - 1;
    for (int g = 2; g < p; g++) {
        int ok = 1;
        for (int i = 0; i < nf; i++) {
            int q = factors_pm1[i].p;
            // Compute g^(m/q) mod p
            long long base = g, exp = m / q, result = 1;
            int mod = p;
            while (exp > 0) {
                if (exp & 1) result = (result * base) % mod;
                base = (base * base) % mod;
                exp >>= 1;
            }
            if (result == 1) { ok = 0; break; }
        }
        if (ok) return g;
    }
    return -1; // should not happen
}

// Discrete log table: dlog[a] = k where root^k == a (mod p), for a=1..p-1.
static int *dlog_table;

static void build_dlog_table(int p, int root) {
    dlog_table = malloc(p * sizeof(int));
    for (int i = 0; i < p; i++) dlog_table[i] = -1;
    long long x = 1;
    int m = p - 1;
    for (int k = 0; k < m; k++) {
        dlog_table[x] = k;
        x = (x * root) % p;
    }
}

// S(p) computation: returns (exact_value, log10_value).
// We store exact_value as mpz and log10 as mpfr.
typedef struct {
    mpz_t value;
    mpfr_t log;
} SResult;

static SResult *s_cache;

static void S_for_prime(int p, mpz_t out_value, mpfr_t out_log) {
    if (s_cache[p].value->_mp_alloc > 0) {
        // Already computed.
        mpz_set(out_value, s_cache[p].value);
        mpfr_set(out_log, s_cache[p].log, MPFR_RNDN);
        return;
    }

    if (p == 2) {
        mpz_set_ui(s_cache[p].value, 1);
        mpfr_set_d(s_cache[p].log, 0.0, MPFR_RNDN);
        mpz_set(out_value, s_cache[p].value);
        mpfr_set(out_log, s_cache[p].log, MPFR_RNDN);
        return;
    }

    int m = p - 1;
    Factor factors_m[32];
    int nf = factorize(m, factors_m);
    int divs[1024];
    int ndivs = divisors_from_factors(factors_m, nf, divs);

    int root = primitive_root(p, factors_m, nf);
    build_dlog_table(p, root);

    // For each proper divisor c of m, find the least rational prime q
    // with gcd(dlog[q % p], m) == c.
    int needed = ndivs - 1; // all divisors except m itself
    int *least_prime_for_c = malloc(ndivs * sizeof(int));
    for (int i = 0; i < ndivs; i++) least_prime_for_c[i] = -1;

    // Map divisor value to index.
    // divs is sorted, so we can binary search.
    int found = 0;
    for (int qi = 0; qi < num_primes && found < needed; qi++) {
        int q = primes[qi];
        if (q == p) continue;
        int r = q % p;
        int dl = dlog_table[r];
        int c = gcd_int(dl, m);
        if (c >= m) continue;
        // Find index of c in divs.
        // Binary search.
        int lo = 0, hi = ndivs - 1;
        int idx = -1;
        while (lo <= hi) {
            int mid = (lo + hi) / 2;
            if (divs[mid] == c) { idx = mid; break; }
            if (divs[mid] < c) lo = mid + 1;
            else hi = mid - 1;
        }
        if (idx >= 0 && least_prime_for_c[idx] < 0) {
            least_prime_for_c[idx] = q;
            found++;
        }
    }

    // best_by_M: for each M in divs, best[d] = least prime q with gcd(c, M) = d.
    // We store as arrays indexed by divisor index.
    // For simplicity, use a 2D structure: best_by_M[M_idx][d_idx] = q.
    // Since arrays are small, use a simple approach.

    // dp_value[h] and dp_log[h], indexed by divisor index.
    mpz_t *dp_value = calloc(ndivs, sizeof(mpz_t));
    mpfr_t *dp_log = calloc(ndivs, sizeof(mpfr_t));
    int *dp_set = calloc(ndivs, sizeof(int));

    for (int i = 0; i < ndivs; i++) {
        mpz_init(dp_value[i]);
        mpfr_init2(dp_log[i], 200);
    }

    // h=1 is at index 0 (divs[0]=1).
    mpz_set_ui(dp_value[0], 1);
    mpfr_set_d(dp_log[0], 0.0, MPFR_RNDN);
    dp_set[0] = 1;

    // Precompute best_by_M for each M.
    // For each M (by index), and each c (by index), compute d = gcd(divs[c_idx], divs[M_idx]).
    // If d < M, store least_prime_for_c[c_idx] if it's the smallest for that d.
    int **best_q = calloc(ndivs, sizeof(int*));
    for (int Mi = 0; Mi < ndivs; Mi++) {
        best_q[Mi] = calloc(ndivs, sizeof(int));
        for (int di = 0; di < ndivs; di++) best_q[Mi][di] = -1;
        int M = divs[Mi];
        if (M == 1) continue;
        for (int ci = 0; ci < ndivs; ci++) {
            if (least_prime_for_c[ci] < 0) continue;
            int c = divs[ci];
            int d = gcd_int(c, M);
            if (d >= M) continue;
            // Find index of d.
            int lo = 0, hi = ndivs - 1, di = -1;
            while (lo <= hi) {
                int mid = (lo + hi) / 2;
                if (divs[mid] == d) { di = mid; break; }
                if (divs[mid] < d) lo = mid + 1;
                else hi = mid - 1;
            }
            if (di >= 0) {
                if (best_q[Mi][di] < 0 || least_prime_for_c[ci] < best_q[Mi][di])
                    best_q[Mi][di] = least_prime_for_c[ci];
            }
        }
    }

    // DP over divisors h.
    for (int hi = 0; hi < ndivs; hi++) {
        if (!dp_set[hi]) continue;
        int h = divs[hi];
        int M = m / h;
        if (M == 1) continue;
        // Find M index.
        int lo = 0, hi2 = ndivs - 1, Mi = -1;
        while (lo <= hi2) {
            int mid = (lo + hi2) / 2;
            if (divs[mid] == M) { Mi = mid; break; }
            if (divs[mid] < M) lo = mid + 1;
            else hi2 = mid - 1;
        }
        if (Mi < 0) continue;

        for (int Li = 0; Li < ndivs; Li++) {
            int L = divs[Li];
            if (L <= 1 || M % L != 0) continue;
            int next_h = h * L;
            // Find index of next_h.
            lo = 0; hi2 = ndivs - 1;
            int nhi = -1;
            while (lo <= hi2) {
                int mid = (lo + hi2) / 2;
                if (divs[mid] == next_h) { nhi = mid; break; }
                if (divs[mid] < next_h) lo = mid + 1;
                else hi2 = mid - 1;
            }
            if (nhi < 0) continue;

            int d = M / L;
            // Find index of d.
            lo = 0; hi2 = ndivs - 1;
            int di = -1;
            while (lo <= hi2) {
                int mid = (lo + hi2) / 2;
                if (divs[mid] == d) { di = mid; break; }
                if (divs[mid] < d) lo = mid + 1;
                else hi2 = mid - 1;
            }
            if (di < 0) continue;

            int q = best_q[Mi][di];
            if (q < 0) continue;

            // candidate = base_value * q^(L-1)
            mpz_t candidate;
            mpz_init(candidate);
            mpz_t qp;
            mpz_init(qp);
            mpz_ui_pow_ui(qp, (unsigned long)q, (unsigned long)(L - 1));
            mpz_mul(candidate, dp_value[hi], qp);

            // candidate_log = base_log + (L-1)*log10(q)
            mpfr_t candidate_log, tmp;
            mpfr_init2(candidate_log, 200);
            mpfr_init2(tmp, 200);
            mpfr_set(candidate_log, dp_log[hi], MPFR_RNDN);
            mpfr_set_d(tmp, log10((double)q), MPFR_RNDN);
            mpfr_mul_ui(tmp, tmp, (unsigned long)(L - 1), MPFR_RNDN);
            mpfr_add(candidate_log, candidate_log, tmp, MPFR_RNDN);

            if (!dp_set[nhi] || mpz_cmp(candidate, dp_value[nhi]) < 0) {
                mpz_set(dp_value[nhi], candidate);
                mpfr_set(dp_log[nhi], candidate_log, MPFR_RNDN);
                dp_set[nhi] = 1;
            }

            mpz_clear(candidate);
            mpz_clear(qp);
            mpfr_clear(candidate_log);
            mpfr_clear(tmp);
        }
    }

    // Find m index (should be last).
    int mi = ndivs - 1; // divs is sorted, m is the largest divisor.
    mpz_set(s_cache[p].value, dp_value[mi]);
    mpfr_set(s_cache[p].log, dp_log[mi], MPFR_RNDN);

    mpz_set(out_value, s_cache[p].value);
    mpfr_set(out_log, s_cache[p].log, MPFR_RNDN);

    // Cleanup.
    for (int i = 0; i < ndivs; i++) {
        mpz_clear(dp_value[i]);
        mpfr_clear(dp_log[i]);
    }
    free(dp_value);
    free(dp_log);
    free(dp_set);
    for (int i = 0; i < ndivs; i++) free(best_q[i]);
    free(best_q);
    free(least_prime_for_c);
    free(dlog_table);
}

// Scientific notation from mpz: mantissa.5f then e then exponent (no +).
void p995_print(void) {
    mpfr_set_default_prec(200);

    sieve(PRIME_SEARCH_LIMIT);

    // Initialize S cache for primes < LIMIT.
    s_cache = calloc(LIMIT, sizeof(SResult));
    for (int i = 0; i < LIMIT; i++) {
        mpz_init(s_cache[i].value);
        mpfr_init2(s_cache[i].log, 200);
    }

    // Compute product of S(p) for all primes p < LIMIT.
    mpz_t product, sp;
    mpz_init_set_ui(product, 1);
    mpz_init(sp);
    mpfr_t sp_log;

    for (int i = 0; i < num_primes; i++) {
        int p = primes[i];
        if (p >= LIMIT) break;
        mpfr_init2(sp_log, 200);
        S_for_prime(p, sp, sp_log);
        mpz_mul(product, product, sp);
        mpfr_clear(sp_log);
    }

    // Convert to scientific notation with 5 decimal places.
    // Get the number of digits.
    size_t num_digits = mpz_sizeinbase(product, 10);
    int exponent = (int)num_digits - 1;

    // Get the first 7 significant digits (6 + 1 for rounding).
    // We need digits[0..6].
    // Extract by dividing product by 10^(num_digits - 7).
    mpz_t divisor, head_z, remainder;
    mpz_init(divisor);
    mpz_init(head_z);
    mpz_init(remainder);

    int significant = 7; // 6 digits + 1 for rounding check
    if ((int)num_digits > significant) {
        mpz_ui_pow_ui(divisor, 10, (unsigned long)(num_digits - significant));
        mpz_tdiv_q(head_z, product, divisor);
    } else {
        mpz_set(head_z, product);
        // Pad: head_z * 10^(significant - num_digits)
        mpz_ui_pow_ui(divisor, 10, (unsigned long)(significant - num_digits));
        mpz_mul(head_z, head_z, divisor);
    }

    // Check rounding: look at the next digit.
    if ((int)num_digits > significant) {
        mpz_tdiv_r(remainder, product, divisor);
        // remainder < divisor, check if remainder * 2 >= divisor
        mpz_mul_ui(remainder, remainder, 2);
        if (mpz_cmp(remainder, divisor) >= 0) {
            mpz_add_ui(head_z, head_z, 1);
            // Check carry
            unsigned long power10 = 1;
            for (int i = 0; i < significant; i++) power10 *= 10;
            if (mpz_cmp_ui(head_z, power10) == 0) {
                mpz_tdiv_q_ui(head_z, head_z, 10);
                exponent++;
            }
        }
    }

    // Format: head_z has 7 digits. mantissa = d.ddddd
    unsigned long head = mpz_get_ui(head_z);
    char buf[64];
    snprintf(buf, sizeof(buf), "%lu", head);
    // buf should have 7 digits.
    int len = strlen(buf);
    // Ensure 7 digits (pad if needed, though shouldn't happen).
    while (len < 7) {
        memmove(buf + 1, buf, len + 1);
        buf[0] = '0';
        len++;
    }

    // mantissa: buf[0].buf[1..5]
    printf("%c.%c%c%c%c%ce%d\n", buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], exponent);

    mpz_clear(product);
    mpz_clear(sp);
    mpz_clear(divisor);
    mpz_clear(head_z);
    mpz_clear(remainder);
}
