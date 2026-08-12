/*
 * Project Euler 763 - Amoebas in a 3D Grid
 *
 * Computes D(10000) mod 1e9 using a recurrence framework based on
 * state compression of the boundary "red path" encoding.
 *
 * The arrays u[n] and v[n] are flat arrays of size n * lens[n],
 * where element (k, idx) is at index k * lens[n] + idx.
 *
 * Ported from the Python reference solver.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define MOD 1000000000LL
#define M_VAL 9999

long long p763_native(void) {
    /* Determine N: largest n with offset[n] <= M, plus 3 */
    int n_tmp = 0;
    while ((n_tmp + 1) * (n_tmp + 2) / 2 <= M_VAL) n_tmp++;
    int max_n = n_tmp - 1;
    int N = max_n + 3;
    int NMAX = N + 2;

    /* offset and lens arrays */
    int *offset_arr = calloc(NMAX, sizeof(int));
    int *lens_arr   = calloc(NMAX, sizeof(int));
    for (int n = 0; n < NMAX; n++) {
        offset_arr[n] = (n + 1) * (n + 2) / 2;
        int ln = M_VAL - offset_arr[n] + 1;
        lens_arr[n] = ln > 0 ? ln : 0;
    }

    /* Allocate u and v arrays (indexed 0..NMAX-1; u[0] unused) */
    uint32_t **u = calloc(NMAX, sizeof(uint32_t *));
    uint32_t **v = calloc(NMAX, sizeof(uint32_t *));
    for (int n = 1; n < NMAX; n++) {
        int ln = lens_arr[n];
        if (ln > 0) {
            u[n] = calloc((size_t)n * ln, sizeof(uint32_t));
            v[n] = calloc((size_t)n * ln, sizeof(uint32_t));
        }
    }

    /* f0 and a2 arrays */
    uint32_t *f0 = calloc(M_VAL + 1, sizeof(uint32_t));
    uint32_t *a2 = calloc(M_VAL + 1, sizeof(uint32_t));
    a2[0] = 1;

    int n_active = 0;
    for (int m = 0; m <= M_VAL; m++) {
        while (n_active + 1 < N + 1 && offset_arr[n_active + 1] <= m)
            n_active++;

        for (int n = 1; n <= n_active; n++) {
            int off = offset_arr[n];
            int ln  = lens_arr[n];
            int idx_cur = m - off;

            int mp1 = m - n - 2;
            int idx1 = mp1 - off;

            int mp2 = m - n - 3;
            int idx2 = mp2 - offset_arr[n + 1];
            int lnp  = lens_arr[n + 1];

            int mp3 = m - n - 1;
            int idx3 = mp3 - offset_arr[n - 1];
            int lnm  = lens_arr[n - 1];

            if (n == 1) {
                int64_t val_u = 0;
                if (idx1 >= 0)
                    val_u += 2LL * u[1][idx1] + v[1][idx1];
                if (idx2 >= 0 && lnp > 0)
                    val_u += (int64_t)v[2][idx2] + u[2][lnp + idx2];
                if (mp3 >= 0)
                    val_u += f0[mp3];
                u[1][idx_cur] = (uint32_t)(val_u % MOD);

                int64_t val_v = 0;
                if (idx1 >= 0)
                    val_v += 2LL * v[1][idx1] + 2LL * u[1][idx1];
                if (idx2 >= 0 && lnp > 0)
                    val_v += (int64_t)v[2][lnp + idx2] + 2LL * u[2][idx2];
                if (mp3 >= 0)
                    val_v += f0[mp3];
                v[1][idx_cur] = (uint32_t)(val_v % MOD);
                continue;
            }

            uint32_t *u_n = u[n], *v_n = v[n];
            uint32_t *u_p = u[n + 1], *v_p = v[n + 1];
            uint32_t *u_m = u[n - 1], *v_m = v[n - 1];

            uint32_t u_n1 = (idx1 >= 0) ? u_n[idx1] : 0;
            uint32_t v_n1 = (idx1 >= 0) ? v_n[idx1] : 0;
            uint32_t u_p1 = (idx2 >= 0 && lnp > 0) ? u_p[idx2] : 0;
            uint32_t v_p1 = (idx2 >= 0 && lnp > 0) ? v_p[idx2] : 0;

            int base = 0, base_next = ln;
            int base_p = lnp, base_m = 0;

            if (idx1 < 0) {
                /* Only n-1 term survives */
                for (int k = 1; k < n; k++) {
                    u_n[base + idx_cur] = u_m[base_m + idx3];
                    v_n[base + idx_cur] = v_m[base_m + idx3];
                    base = base_next;
                    base_next += ln;
                    base_m += lnm;
                }
                u_n[(n - 1) * ln + idx_cur] = u_m[(n - 2) * lnm + idx3];
                v_n[(n - 1) * ln + idx_cur] = v_m[(n - 2) * lnm + idx3];
                continue;
            }

            if (idx2 >= 0 && lnp > 0) {
                /* Full recurrence */
                for (int k = 1; k < n; k++) {
                    int64_t uu = (int64_t)u_n[base + idx1]
                        + v_p1
                        + u_p[base_p + idx2]
                        + u_m[base_m + idx3]
                        + v_n1
                        + u_n[base_next + idx1];
                    u_n[base + idx_cur] = (uint32_t)(uu % MOD);

                    int64_t vv = (int64_t)v_n[base + idx1]
                        + v_p[base_p + idx2]
                        + u_p1
                        + v_m[base_m + idx3]
                        + v_n[base_next + idx1]
                        + u_n1;
                    v_n[base + idx_cur] = (uint32_t)(vv % MOD);

                    base = base_next;
                    base_next += ln;
                    base_p += lnp;
                    base_m += lnm;
                }
                int bl = (n - 1) * ln;
                int64_t uu = 2LL * u_n[bl + idx1]
                    + v_n1
                    + v_p1
                    + u_p[base_p + idx2]
                    + u_m[(n - 2) * lnm + idx3];
                u_n[bl + idx_cur] = (uint32_t)(uu % MOD);

                int64_t vv = 2LL * v_n[bl + idx1]
                    + 2LL * u_n1
                    + v_p[base_p + idx2]
                    + 2LL * u_p1
                    + v_m[(n - 2) * lnm + idx3];
                v_n[bl + idx_cur] = (uint32_t)(vv % MOD);
            } else {
                /* No n+1 term */
                for (int k = 1; k < n; k++) {
                    int64_t uu = (int64_t)u_n[base + idx1]
                        + u_m[base_m + idx3]
                        + v_n1
                        + u_n[base_next + idx1];
                    u_n[base + idx_cur] = (uint32_t)(uu % MOD);

                    int64_t vv = (int64_t)v_n[base + idx1]
                        + v_m[base_m + idx3]
                        + v_n[base_next + idx1]
                        + u_n1;
                    v_n[base + idx_cur] = (uint32_t)(vv % MOD);

                    base = base_next;
                    base_next += ln;
                    base_m += lnm;
                }
                int bl = (n - 1) * ln;
                int64_t uu = 2LL * u_n[bl + idx1]
                    + v_n1
                    + u_m[(n - 2) * lnm + idx3];
                u_n[bl + idx_cur] = (uint32_t)(uu % MOD);

                int64_t vv = 2LL * v_n[bl + idx1]
                    + 2LL * u_n1
                    + v_m[(n - 2) * lnm + idx3];
                v_n[bl + idx_cur] = (uint32_t)(vv % MOD);
            }
        }

        /* f0 and a2 */
        int64_t val_f = 0;
        if (m - 1 >= 0) val_f += a2[m - 1];
        if (m - 2 >= 0) val_f += 4LL * f0[m - 2];
        int mp = m - 3;
        if (mp >= offset_arr[1] && lens_arr[1] > 0) {
            int id1 = mp - offset_arr[1];
            val_f += 2LL * u[1][id1] + v[1][id1];
        }
        f0[m] = (uint32_t)(val_f % MOD);

        if (m >= 1) {
            int64_t val_a = 3LL * a2[m - 1];
            if (m - 2 >= 0) val_a += 3LL * f0[m - 2];
            a2[m] = (uint32_t)(val_a % MOD);
        }
    }

    int64_t ans = a2[M_VAL];

    /* Cleanup */
    for (int n = 1; n < NMAX; n++) {
        free(u[n]);
        free(v[n]);
    }
    free(u); free(v);
    free(f0); free(a2);
    free(offset_arr); free(lens_arr);

    return ans;
}
