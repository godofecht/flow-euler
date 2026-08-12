#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;

static const i64 MOD = 1234567891LL;

static i64 mod_pow(i64 a, i64 e, i64 mod) {
    i64 r = 1 % mod; a %= mod;
    if (a < 0) a += mod;
    while (e > 0) {
        if (e & 1) r = (i64)((__int128)r * a % mod);
        a = (i64)((__int128)a * a % mod);
        e >>= 1;
    }
    return r;
}

/* Batch-invert consecutive integers start, start+1, ..., start+length-1 mod p. */
static void invert_consecutive(i64 start, i64 length, i64 mod, i64 *out) {
    if (length <= 0) return;
    i64 *pref = (i64 *)malloc(length * sizeof(i64));
    i64 acc = 1;
    for (i64 i = 0; i < length; i++) {
        acc = (i64)((__int128)acc * ((start + i) % mod) % mod);
        pref[i] = acc;
    }
    i64 inv_acc = mod_pow(pref[length - 1], mod - 2, mod);
    for (i64 i = length - 1; i >= 0; i--) {
        i64 prev = (i == 0) ? 1 : pref[i - 1];
        out[i] = (i64)((__int128)inv_acc * prev % mod);
        inv_acc = (i64)((__int128)inv_acc * ((start + i) % mod) % mod);
    }
    free(pref);
}

/* Batch-invert a list of values mod p. */
static void invert_list(i64 *vals, i64 n, i64 mod, i64 *out) {
    if (n == 0) return;
    i64 *pref = (i64 *)malloc(n * sizeof(i64));
    i64 acc = 1;
    for (i64 i = 0; i < n; i++) {
        acc = (i64)((__int128)acc * vals[i] % mod);
        pref[i] = acc;
    }
    i64 inv_acc = mod_pow(pref[n - 1], mod - 2, mod);
    for (i64 i = n - 1; i >= 0; i--) {
        i64 prev = (i == 0) ? 1 : pref[i - 1];
        out[i] = (i64)((__int128)inv_acc * prev % mod);
        inv_acc = (i64)((__int128)inv_acc * vals[i] % mod);
    }
    free(pref);
}

/* C(n,k) mod prime, product formula with batch inversion. */
static i64 binom_mod(i64 n, i64 k, i64 mod) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    if (k == 0) return 1;
    i64 base = n - k;
    i64 res = 1;
    i64 start = 1;
    i64 block = 200000;
    i64 *invs = (i64 *)malloc(block * sizeof(i64));
    while (start <= k) {
        i64 length = block;
        if (length > k - start + 1) length = k - start + 1;
        invert_consecutive(start, length, mod, invs);
        i64 b = base + start;
        for (i64 i = 0; i < length; i++) {
            res = (i64)((__int128)res * ((b + i) % mod) % mod);
            res = (i64)((__int128)res * invs[i] % mod);
        }
        start += length;
    }
    free(invs);
    return res;
}

static i64 amidakuji_count_mod(i64 m, i64 n, i64 mod) {
    i64 L = m + n;
    if (L & 1) return 0;

    i64 t = L / 2;
    i64 k = m & 1;
    i64 limit = (m < n) ? m : n;

    i64 layout;
    if (k == 0) {
        layout = binom_mod(t, m / 2, mod);
    } else {
        layout = (i64)((__int128)(t % mod) * binom_mod(t - 1, (m - 1) / 2, mod) % mod);
    }

    i64 inv3 = mod_pow(3, mod - 2, mod);
    i64 pow2 = (k == 0) ? 1 : 2;
    i64 sign = k ? (mod - 1) : 1;

    i64 total = 0;
    i64 block = 200000;
    i64 *nums = (i64 *)malloc(block * sizeof(i64));
    i64 *dens = (i64 *)malloc(block * sizeof(i64));
    i64 *inv_dens = (i64 *)malloc(block * sizeof(i64));

    while (k <= limit) {
        i64 steps = block;
        i64 avail = ((limit - k) / 2) + 1;
        if (steps > avail) steps = avail;

        i64 kk = k;
        for (i64 i = 0; i < steps; i++) {
            nums[i] = (i64)((__int128)((m - kk) % mod) * ((n - kk) % mod) % mod);
            dens[i] = (i64)((__int128)((4 * (kk + 1)) % mod) * ((kk + 2) % mod) % mod);
            kk += 2;
        }

        invert_list(dens, steps, mod, inv_dens);

        for (i64 i = 0; i < steps; i++) {
            i64 orientations = (i64)((__int128)((pow2 + 2 * sign) % mod) * inv3 % mod);
            total = (i64)((__int128)(total + (i64)((__int128)layout * orientations % mod)) % mod);
            layout = (i64)((__int128)((i64)((__int128)layout * nums[i] % mod)) * inv_dens[i] % mod);
            pow2 = (i64)((__int128)pow2 * 4 % mod);
        }

        k += 2 * steps;
    }

    free(nums);
    free(dens);
    free(inv_dens);
    return total;
}

long long p837_native(void) {
    return (long long)amidakuji_count_mod(123456789LL, 987654321LL, MOD);
}
