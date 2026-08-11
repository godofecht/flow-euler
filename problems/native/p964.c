/*
 * Project Euler 964
 *
 * Computes P(k) for k = 7 using exact rational arithmetic (GMP).
 *
 *   n = k*(k-1)/2 + 1
 *   P(k) = (1/n!) * sum_{m=0}^{n-1} (-1)^m * prod_{i=1..k} mult_{lam,i}
 *          / dim(lam)^(k-1)
 *
 * where lam is the hook partition (n-m, 1^m), dim = C(n-1, m), and
 * mult_{lam,i} counts corner-removal paths from lam down to (i).
 *
 * Reference: /tmp/pes_ref/solvers/964.py
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --- memoized corner-removal path counts ------------------------------- */

#define TBITS (1 << 19)
#define TMASK (TBITS - 1)

typedef struct {
    char *key;
    mpz_t val;
} Entry;

static Entry table[TBITS];

static size_t hash_str(const char *s) {
    size_t h = 1469598103934665603ULL;
    while (*s) {
        h ^= (unsigned char)*s;
        h *= 1099511628211ULL;
        s++;
    }
    return h;
}

/* Encode partition (parts >= 1, so no NUL bytes) plus target i. */
static void make_key(char *buf, const int *start, int len, int i) {
    buf[0] = (char)i;
    for (int j = 0; j < len; j++) buf[1 + j] = (char)start[j];
    buf[1 + len] = '\0';
}

static int lookup(const char *key, mpz_t out) {
    size_t h = hash_str(key) & TMASK;
    while (table[h].key) {
        if (strcmp(table[h].key, key) == 0) {
            mpz_set(out, table[h].val);
            return 1;
        }
        h = (h + 1) & TMASK;
    }
    return 0;
}

static void insert(const char *key, const mpz_t val) {
    size_t h = hash_str(key) & TMASK;
    while (table[h].key) {
        if (strcmp(table[h].key, key) == 0) {
            mpz_set(table[h].val, val);
            return;
        }
        h = (h + 1) & TMASK;
    }
    table[h].key = strdup(key);
    mpz_init_set(table[h].val, val);
}

/* Count paths removing corner boxes from (start, len) until reaching (i). */
static void count_paths(const int *start, int len, int i, mpz_t out) {
    int sum = 0;
    for (int j = 0; j < len; j++) sum += start[j];
    if (sum < i) { mpz_set_ui(out, 0); return; }
    if (len == 1 && start[0] == i) { mpz_set_ui(out, 1); return; }

    char key[64];
    make_key(key, start, len, i);
    if (lookup(key, out)) return;

    mpz_t total, tmp;
    mpz_init_set_ui(total, 0);
    mpz_init(tmp);

    int newp[32];
    for (int r = 0; r < len; r++) {
        if (r == len - 1 || start[r + 1] < start[r]) {
            int nlen = len;
            for (int j = 0; j < len; j++) newp[j] = start[j];
            newp[r] -= 1;
            if (newp[r] == 0) {
                for (int j = r; j < nlen - 1; j++) newp[j] = newp[j + 1];
                nlen--;
            }
            if (nlen == 0) continue;
            if (newp[0] < i) continue;       /* must still contain (i) */
            count_paths(newp, nlen, i, tmp);
            mpz_add(total, total, tmp);
        }
    }

    mpz_set(out, total);
    insert(key, out);
    mpz_clear(total);
    mpz_clear(tmp);
}

/* --- main computation -------------------------------------------------- */

double p964_native(void) {
    const int k = 7;
    const int n = k * (k - 1) / 2 + 1;       /* 22 */

    mpz_t nfact;
    mpz_init(nfact);
    mpz_fac_ui(nfact, n);                    /* n! */

    mpq_t total, term;
    mpq_init(total);
    mpq_init(term);

    mpz_t prod, signed_prod, d, d6, mult;
    mpz_init(prod);
    mpz_init(signed_prod);
    mpz_init(d);
    mpz_init(d6);
    mpz_init(mult);

    for (int m = 0; m < n; m++) {
        /* hook partition (n-m, 1^m) */
        int lam[32];
        int lamlen;
        if (m == 0) {
            lam[0] = n;
            lamlen = 1;
        } else {
            lam[0] = n - m;
            for (int j = 1; j <= m; j++) lam[j] = 1;
            lamlen = m + 1;
        }

        /* prod = product of multiplicities for i = 1..k */
        mpz_set_ui(prod, 1);
        for (int i = 1; i <= k; i++) {
            count_paths(lam, lamlen, i, mult);
            mpz_mul(prod, prod, mult);
            if (mpz_sgn(prod) == 0) break;
        }

        /* dim = C(n-1, m) */
        mpz_bin_uiui(d, n - 1, m);

        /* d6 = dim^(k-1) = dim^6 */
        mpz_pow_ui(d6, d, k - 1);

        /* term = (-1)^m * prod / d6 */
        mpz_set(signed_prod, prod);
        if (m % 2) mpz_neg(signed_prod, signed_prod);
        mpz_set(mpq_numref(term), signed_prod);
        mpz_set(mpq_denref(term), d6);
        mpq_canonicalize(term);
        mpq_add(total, total, term);
    }

    /* P = total / n! */
    mpz_t pden;
    mpz_init(pden);
    mpz_mul(pden, mpq_denref(total), nfact);
    mpq_t P;
    mpq_init(P);
    mpz_set(mpq_numref(P), mpq_numref(total));
    mpz_set(mpq_denref(P), pden);
    mpq_canonicalize(P);

    double res = mpq_get_d(P);

    mpz_clear(nfact);
    mpq_clear(total);
    mpq_clear(term);
    mpz_clear(prod);
    mpz_clear(signed_prod);
    mpz_clear(d);
    mpz_clear(d6);
    mpz_clear(mult);
    mpz_clear(pden);
    mpq_clear(P);

    return res;
}
