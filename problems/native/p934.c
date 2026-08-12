/* Project Euler 934: Unlucky Primes
 * U(N) = sum_{n=1..N} u(n) where u(n) is the smallest prime p
 * such that (n mod p) is NOT a multiple of 7.
 * Compute U(10^17) using CRT residue construction + explicit enumeration.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef int64_t i64;
typedef uint64_t u64;
typedef __int128 i128;

/* Extended Euclid for modular inverse */
static i64 inv_mod(i64 a, i64 m) {
    a %= m;
    if (a < 0) a += m;
    i64 t0 = 0, t1 = 1;
    i64 r0 = m, r1 = a;
    while (r1) {
        i64 q = r0 / r1;
        i64 tmp;
        tmp = t0 - q * t1; t0 = t1; t1 = tmp;
        tmp = r0 - q * r1; r0 = r1; r1 = tmp;
    }
    if (r0 != 1) return -1; /* should not happen */
    t0 %= m;
    if (t0 < 0) t0 += m;
    return t0;
}

/* Simple prime sieve */
static int is_prime(int x) {
    if (x < 2) return 0;
    if (x == 2) return 1;
    if (x % 2 == 0) return 0;
    for (int i = 3; (i64)i * i <= x; i += 2) {
        if (x % i == 0) return 0;
    }
    return 1;
}

long long p934_native(void) {
    i64 N = 1;
    for (int i = 0; i < 17; i++) N *= 10; /* 10^17 */

    i64 ans = 0;
    i64 c_prev = N;

    /* Phase 1: CRT mode */
    i64 M = 1;
    i64 *residues = (i64*)malloc(sizeof(i64) * 1);
    residues[0] = 0;
    int num_residues = 1;
    int cap_residues = 1;

    /* Phase 2: explicit survivors */
    i64 *survivors = NULL;
    int num_survivors = 0;
    int cap_survivors = 0;

    for (int p = 2; ; p++) {
        if (!is_prime(p)) continue;

        /* Allowed residues mod p: multiples of 7 in [0, p-1] */
        int num_allowed = (p + 6) / 7; /* ceil(p/7) = number of multiples of 7 in [0, p-1] */
        /* Actually range(0, p, 7) gives 0, 7, 14, ... < p. Count = (p-1)/7 + 1 = (p+6)/7 */
        /* But for p=2: (2+6)/7 = 1 (just 0). For p=7: (7+6)/7 = 1 (just 0). For p=11: (11+6)/7 = 2 (0, 7). Correct. */

        i64 c;

        if (survivors == NULL) {
            /* CRT mode */
            i128 M_new = (i128)M * p;

            if (M_new <= (i128)N) {
                /* Extend residues using CRT */
                i64 inv = inv_mod(M, p);
                i64 M_old = M;

                int new_cap = num_residues * num_allowed;
                i64 *new_residues = (i64*)malloc(sizeof(i64) * new_cap);
                int nr = 0;

                for (int i = 0; i < num_residues; i++) {
                    i64 r = residues[i];
                    for (int s = 0; s < p; s += 7) {
                        i64 diff = s - r;
                        i64 t = ((diff % p) * inv) % p;
                        if (t < 0) t += p;
                        i64 new_r = r + M_old * t;
                        new_residues[nr++] = new_r;
                    }
                }

                free(residues);
                residues = new_residues;
                num_residues = nr;
                cap_residues = new_cap;
                M = (i64)M_new;

                /* Count numbers <= N matching these residues (period M) */
                i64 q = N / M;
                i64 rem = N % M;
                i64 extra = 0;
                for (int i = 0; i < num_residues; i++) {
                    if (residues[i] > 0 && residues[i] <= rem) {
                        extra++;
                    }
                }
                c = q * num_residues + extra;

            } else {
                /* Switch to explicit enumeration */
                i64 inv = inv_mod(M, p);
                i64 M_old = M;

                cap_survivors = num_residues * num_allowed;
                survivors = (i64*)malloc(sizeof(i64) * cap_survivors);
                num_survivors = 0;

                for (int i = 0; i < num_residues; i++) {
                    i64 r = residues[i];
                    for (int s = 0; s < p; s += 7) {
                        i64 diff = s - r;
                        i64 t = ((diff % p) * inv) % p;
                        if (t < 0) t += p;
                        i64 x = r + M_old * t;
                        if (x > 0 && x <= N) {
                            survivors[num_survivors++] = x;
                        }
                    }
                }

                free(residues);
                residues = NULL;
                M = (i64)M_new;
                c = num_survivors;
            }

        } else {
            /* Explicit enumeration mode: filter survivors */
            int j = 0;
            for (int i = 0; i < num_survivors; i++) {
                i64 n = survivors[i];
                if ((n % p) % 7 == 0) {
                    survivors[j++] = n;
                }
            }
            num_survivors = j;
            c = num_survivors;
        }

        /* Numbers that stop at prime p have u(n) = p */
        ans += (i64)p * (c_prev - c);
        c_prev = c;

        if (c == 0) break;
    }

    if (residues) free(residues);
    if (survivors) free(survivors);

    return ans;
}
