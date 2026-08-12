// Project Euler 902: Permutation Powers
// Compute P(100) mod 1e9+7.
// Port of the Python reference solver.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MOD 1000000007LL
#define A 1000000007LL  // multiplier in tau
#define MM 100          // m = 100

static long long mod_pow(long long base, long long exp, long long mod) {
    long long result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = result * base % mod;
        base = base * base % mod;
        exp >>= 1;
    }
    return result;
}

static long long egcd_inv(long long a, long long b) {
    long long x0 = 1, y0 = 0, x1 = 0, y1 = 1;
    while (b) {
        long long q = a / b;
        long long r = a - q * b;
        a = b; b = r;
        long long tx = x0 - q * x1;
        x0 = x1; x1 = tx;
    }
    return ((x0 % MOD) + MOD) % MOD;  // not used directly, see below
}

// Extended GCD to find modular inverse of a mod m
static long long inv_mod(long long a, long long m) {
    long long x0 = 1, y0 = 0, x1 = 0, y1 = 1;
    long long aa = a % m;
    long long bb = m;
    while (bb) {
        long long q = aa / bb;
        long long r = aa - q * bb;
        aa = bb; bb = r;
        long long tx = x0 - q * x1;
        x0 = x1; x1 = tx;
    }
    // aa = gcd, x0 = inverse
    return ((x0 % m) + m) % m;
}

static long long gcd_ll(long long a, long long b) {
    while (b) {
        long long t = a % b;
        a = b; b = t;
    }
    return a;
}

// Compute L = lcm(1..m) mod MOD
static long long lcm_upto_mod(int m) {
    // L = product of p^floor(log_p(m)) for primes p <= m
    long long result = 1;
    char *is_prime = (char *)malloc(m + 1);
    memset(is_prime, 1, m + 1);
    is_prime[0] = is_prime[1] = 0;
    for (int i = 2; i <= m; i++) {
        if (is_prime[i]) {
            // Find largest power of i <= m
            long long pk = 1;
            while (pk * i <= m) pk *= i;
            result = result * (pk % MOD) % MOD;
            for (int j = i * i; j <= m; j += i) is_prime[j] = 0;
        }
    }
    free(is_prime);
    return result;
}

// Binary search: lower_bound (first index where arr[idx] >= val)
static int lower_bound(int *arr, int len, int val) {
    int lo = 0, hi = len;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (arr[mid] < val) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

static int cmp_int(const void *a, const void *b) {
    return *(const int *)a - *(const int *)b;
}

long long p902_native(void) {
    int m = MM;
    int n = m * (m + 1) / 2;  // 5050

    // Precompute factorials mod MOD
    long long *fact_n = (long long *)malloc((n + 1) * sizeof(long long));
    fact_n[0] = fact_n[1] = 1;
    for (int i = 2; i <= n; i++) {
        fact_n[i] = fact_n[i - 1] * i % MOD;
    }

    // weights[i] = fact_n[n - i] for i = 1..n
    long long *weights = (long long *)malloc((n + 1) * sizeof(long long));
    for (int i = 1; i <= n; i++) {
        weights[i] = fact_n[n - i];
    }

    // Build cycles
    // cycles[l] = array of l ints
    int **cycles = (int **)calloc(m + 1, sizeof(int *));
    int *clen = (int *)calloc(n + 1, sizeof(int));
    int *coff = (int *)calloc(n + 1, sizeof(int));

    long long inva = inv_mod(A % n, n);

    for (int l = 1; l <= m; l++) {
        int start = l * (l - 1) / 2 + 1;
        int *cyc = (int *)malloc(l * sizeof(int));
        for (int x = start; x < start + l; x++) {
            // tau_inv(x): t = (inva * ((x-1) % n)) % n; return n if t==0 else t
            long long t = (inva * ((long long)(x - 1) % n)) % n;
            cyc[x - start] = (t == 0) ? n : (int)t;
        }
        cycles[l] = cyc;
        for (int idx = 0; idx < l; idx++) {
            int elem = cyc[idx];
            clen[elem] = l;
            coff[elem] = idx;
        }
    }

    // L mod MOD
    long long L_mod = lcm_upto_mod(m);

    // Precompute comp, gcd_tab, scale
    // comp[a][b][d] for a,b in 1..m, d in 0..gcd(a,b)-1
    // Use flat arrays
    long long gcd_tab[MM + 1][MM + 1];
    long long scale[MM + 1][MM + 1];
    // comp stored as comp[a][b] = pointer to array of g ints
    int *comp[MM + 1][MM + 1];
    memset(comp, 0, sizeof(comp));

    for (int a = 1; a <= m; a++) {
        for (int b = 1; b <= m; b++) {
            int g = (int)gcd_ll(a, b);
            gcd_tab[a][b] = g;
            int l = (a / g) * b;  // lcm(a,b)
            scale[a][b] = L_mod * mod_pow(l % MOD, MOD - 2, MOD) % MOD;

            int *Avals = cycles[a];
            int *Bvals = cycles[b];

            // Partition Avals by residue mod g: A_subs[r] = elements at indices r, r+g, r+2g, ...
            // B_sorted[r] = sorted elements at indices r, r+g, ...
            int a_len = a, b_len = b;
            int a_sub_count[MM];  // number of elements in each residue class
            int b_sub_count[MM];
            for (int r = 0; r < g; r++) {
                a_sub_count[r] = (a_len - r + g - 1) / g;  // ceil((a_len - r) / g)
                if (r >= a_len) a_sub_count[r] = 0;
                b_sub_count[r] = (b_len - r + g - 1) / g;
                if (r >= b_len) b_sub_count[r] = 0;
            }

            // Build A_subs and B_sorted
            int *A_subs[MM];
            int *B_sorted[MM];
            for (int r = 0; r < g; r++) {
                A_subs[r] = (int *)malloc(a_sub_count[r] * sizeof(int) + 1);
                B_sorted[r] = (int *)malloc(b_sub_count[r] * sizeof(int) + 1);
                int cnt = 0;
                for (int idx = r; idx < a_len; idx += g) {
                    A_subs[r][cnt++] = Avals[idx];
                }
                cnt = 0;
                for (int idx = r; idx < b_len; idx += g) {
                    B_sorted[r][cnt++] = Bvals[idx];
                }
                qsort(B_sorted[r], b_sub_count[r], sizeof(int), cmp_int);
            }

            int *smallcounts = (int *)calloc(g, sizeof(int));
            for (int d = 0; d < g; d++) {
                int tot = 0;
                for (int r = 0; r < g; r++) {
                    int s = ((r - d) % g + g) % g;
                    int *Bs = B_sorted[s];
                    int bs_len = b_sub_count[s];
                    for (int idx2 = 0; idx2 < a_sub_count[r]; idx2++) {
                        int aval = A_subs[r][idx2];
                        tot += lower_bound(Bs, bs_len, aval);
                    }
                }
                smallcounts[d] = tot;
            }
            comp[a][b] = smallcounts;

            for (int r = 0; r < g; r++) {
                free(A_subs[r]);
                free(B_sorted[r]);
            }
        }
    }

    // Sum of ranks over all exponents k mod L
    long long acc = 0;
    for (int i = 1; i < n; i++) {
        long long wi = weights[i];
        int ai = clen[i];
        int offi = coff[i];
        for (int j = i + 1; j <= n; j++) {
            int aj = clen[j];
            int g = (int)gcd_tab[ai][aj];
            int d;
            if (g == 1) {
                d = 0;
            } else {
                d = ((offi - coff[j]) % g + g) % g;
            }
            int small = comp[ai][aj][d];
            if (small) {
                long long cnt = scale[ai][aj] * small % MOD;
                acc += wi * cnt;
                if (acc >= (1LL << 62)) acc %= MOD;
            }
        }
    }
    acc %= MOD;
    long long sum_ranks_mod = (L_mod + acc) % MOD;

    // m! mod MOD
    long long m_fact = 1;
    for (int k = 2; k <= m; k++) {
        m_fact = m_fact * k % MOD;
    }

    long long ans = sum_ranks_mod * m_fact % MOD * mod_pow(L_mod, MOD - 2, MOD) % MOD;

    // Cleanup
    for (int l = 1; l <= m; l++) free(cycles[l]);
    free(cycles); free(clen); free(coff); free(fact_n); free(weights);
    for (int a = 1; a <= m; a++)
        for (int b = 1; b <= m; b++)
            free(comp[a][b]);

    return ans;
}
