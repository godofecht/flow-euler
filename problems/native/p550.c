
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef int64_t i64;

#define MOD 987654321LL

static i64 mod_pow(i64 a, i64 e, i64 mod) {
    i64 r = 1 % mod;
    a %= mod;
    if (a < 0) a += mod;
    while (e > 0) {
        if (e & 1) r = (i64)((__int128)r * a % mod);
        a = (i64)((__int128)a * a % mod);
        e >>= 1;
    }
    return r;
}

/* Extended Euclidean modular inverse (MOD is not prime, so Fermat is invalid). */
static i64 mod_inv(i64 a, i64 mod) {
    i64 old_r = a % mod, r = mod;
    i64 old_s = 1, s = 0;
    if (old_r < 0) old_r += mod;
    while (r != 0) {
        i64 q = old_r / r;
        i64 t = old_r - q * r; old_r = r; r = t;
        i64 t2 = old_s - q * s; old_s = s; s = t2;
    }
    /* old_r is gcd, must be 1 for a coprime to mod. */
    if (old_s < 0) old_s += mod;
    return old_s % mod;
}

/* In-place Walsh-Hadamard transform for XOR convolution. */
static void fwht_xor(i64 *a, i64 n, i64 mod) {
    for (i64 h = 1; h < n; h <<= 1) {
        i64 step = h << 1;
        for (i64 i = 0; i < n; i += step) {
            for (i64 j = i; j < i + h; j++) {
                i64 x = a[j];
                i64 y = a[j + h];
                i64 s = x + y;
                if (s >= mod) s -= mod;
                i64 d = x - y;
                if (d < 0) d += mod;
                a[j] = s;
                a[j + h] = d;
            }
        }
    }
}

/*
 * h[t] = nimber for any pile with Omega(n)=t.
 * h[t] = mex( {h[i] xor h[j] : 1<=i,j<t} ).
 * tmax <= 23 for n=1e7, so this is tiny.
 * Returns max nimber value produced.
 */
static i64 compute_h_sequence(i64 tmax, i64 *h) {
    if (tmax <= 0) return 0;
    /* h[0] unused, h[1] = 0 */
    h[1] = 0;
    i64 max_g = 0;
    for (i64 t = 2; t <= tmax; t++) {
        /* Collect reachable nimbers h[i]^h[j] for 1<=i,j<t. */
        /* tmax <= 23, nimbers are small; use a bitmap. */
        uint64_t reachable[64];
        for (i64 z = 0; z < 64; z++) reachable[z] = 0;
        for (i64 i = 1; i < t; i++) {
            for (i64 j = 1; j < t; j++) {
                i64 v = h[i] ^ h[j];
                reachable[(uint64_t)v >> 6] |= (1ULL << (v & 63));
            }
        }
        i64 mex = 0;
        while (mex < 64 * 64) {
            if (reachable[(uint64_t)mex >> 6] & (1ULL << (mex & 63))) {
                mex++;
            } else {
                break;
            }
        }
        h[t] = mex;
        if (mex > max_g) max_g = mex;
    }
    return max_g;
}

/*
 * Linear sieve: count integers m in [2..n] by Omega(m).
 * counts[t] = how many m in [2..n] have Omega(m)=t.
 * Returns max omega observed.
 */
static i64 omega_counts_up_to(i64 n, i64 *counts) {
    /* spf: smallest prime factor (int32). omega: byte. */
    int32_t *spf = (int32_t *)calloc((size_t)(n + 1), sizeof(int32_t));
    unsigned char *omega = (unsigned char *)calloc((size_t)(n + 1), 1);
    /* ~664579 primes below 1e7; round up. */
    int32_t *primes = (int32_t *)malloc(sizeof(int32_t) * 700000);
    if (!spf || !omega || !primes) {
        free(spf); free(omega); free(primes);
        return -1;
    }
    for (i64 z = 0; z < 32; z++) counts[z] = 0;

    i64 np = 0;
    i64 max_om = 0;

    for (i64 i = 2; i <= n; i++) {
        int32_t si = spf[i];
        i64 oi;
        if (si == 0) {
            spf[i] = (int32_t)i;
            primes[np++] = (int32_t)i;
            omega[i] = 1;
            oi = 1;
            si = (int32_t)i;
        } else {
            oi = (i64)omega[i];
        }
        counts[oi]++;
        if (oi > max_om) max_om = oi;

        for (i64 j = 0; j < np; j++) {
            i64 p = (i64)primes[j];
            i64 ip = i * p;
            if (ip > n) break;
            spf[ip] = (int32_t)p;
            omega[ip] = (unsigned char)(oi + 1);
            if (p == (i64)si) break;
        }
    }

    free(spf);
    free(omega);
    free(primes);
    return max_om;
}

static i64 f(i64 n, i64 k) {
    i64 counts[32];
    i64 max_om = omega_counts_up_to(n, counts);
    if (max_om < 0) return 0;

    i64 h[32];
    for (i64 z = 0; z < 32; z++) h[z] = 0;
    i64 max_g = compute_h_sequence(max_om, h);

    /* size = power of two > max_g */
    i64 size = 1;
    while (size <= max_g) size <<= 1;

    i64 *vec = (i64 *)calloc((size_t)size, sizeof(i64));
    if (!vec) return 0;

    for (i64 t = 1; t <= max_om; t++) {
        vec[h[t]] = (vec[h[t]] + counts[t]) % MOD;
    }

    fwht_xor(vec, size, MOD);

    i64 s = 0;
    for (i64 i = 0; i < size; i++) {
        s = (s + mod_pow(vec[i], k, MOD)) % MOD;
    }

    i64 losing = (i64)((__int128)s * mod_inv(size, MOD) % MOD);
    i64 total = mod_pow(n - 1, k, MOD);
    i64 ans = (total - losing) % MOD;
    if (ans < 0) ans += MOD;

    free(vec);
    return ans;
}

long long p550_native(void) {
    return (long long)f(10000000LL, 1000000000000LL);
}
