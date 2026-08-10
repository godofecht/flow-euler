#include <gmp.h>
#include <stdio.h>
#include <string.h>
typedef struct { int len; int digs[20]; } Part;
static Part parts[10][200];
static int part_count[10];
static void build_partitions(void) {
    memset(part_count, 0, sizeof(part_count));
    parts[0][part_count[0]++].len = 0;
    for (int opt = 1; opt <= 9; opt++) {
        for (int i = 0; i <= 9 - opt; i++) {
            int pc = part_count[i];
            for (int p = 0; p < pc; p++) {
                Part np = parts[i][p];
                /* prepend like python: (options,) + y */
                for (int j = np.len; j > 0; j--) np.digs[j] = np.digs[j-1];
                np.digs[0] = opt;
                np.len++;
                parts[i + opt][part_count[i + opt]++] = np;
            }
        }
    }
}
static void multinomial(mpz_t out, int n, int *d) {
    mpz_set_ui(out, 1);
    int rem = n;
    for (int v = 0; v < 10; v++) {
        for (int i = 1; i <= d[v]; i++) {
            mpz_mul_ui(out, out, (unsigned long)rem);
            mpz_divexact_ui(out, out, (unsigned long)i);
            rem--;
        }
    }
}
long long pe_solve(void) {
    int n = 2020;
    build_partitions();
        mpz_t rep, total, v, comb, mod;
    mpz_inits(rep, total, v, comb, mod, NULL);
    mpz_ui_pow_ui(mod, 10, 16);
    mpz_set_ui(rep, 0);
    for (int i = 0; i < n; i++) { mpz_mul_ui(rep, rep, 10); mpz_add_ui(rep, rep, 1); }
    mpz_set_ui(total, 0);
    for (int k = 0; k <= 9; k++) {
        for (int pi = 0; pi < part_count[k]; pi++) {
            Part p = parts[k][pi];
            p.digs[p.len++] = k;
            if (p.len - 1 >= n) continue;
            int d[10] = {0};
            for (int i = 0; i < p.len; i++) d[p.digs[i]]++;
            int sumd = 0; for (int vv = 1; vv <= 9; vv++) sumd += d[vv];
            d[0] = n - sumd;
            if (d[0] < 0) continue;
            multinomial(comb, n, d);
            mpz_mul(v, rep, comb);
            mpz_mul_ui(v, v, (unsigned long)(2 * k));
            mpz_divexact_ui(v, v, (unsigned long)n);
            mpz_add(total, total, v);
            mpz_mod(total, total, mod);
        }
    }
    long long ans = (long long)mpz_get_ui(total);
    mpz_clears(rep, total, v, comb, mod, NULL);
    return ans;
}
