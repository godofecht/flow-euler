// Project Euler 873: Words with Gaps
// W(10^6, 10^7, 10^8) mod 1e9+7
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
typedef unsigned long long u64;
typedef __int128 i128;

#define MOD 1000000007LL

static i64 modpow(i64 base, i64 exp, i64 mod) {
    i64 r = 1; base %= mod; if (base < 0) base += mod;
    while (exp > 0) {
        if (exp & 1) r = (i64)((i128)r * base % mod);
        base = (i64)((i128)base * base % mod);
        exp >>= 1;
    }
    return r;
}

// Batch modular inverses of [start..end] inclusive.
// Returns malloc'd array of size (end-start+1), invs[i] = (start+i)^(-1) mod mod.
static i64 *batch_inverses(i64 start, i64 end, i64 mod) {
    i64 length = end - start + 1;
    i64 *pref = (i64*)malloc((size_t)(length + 1) * sizeof(i64));
    i64 *invs = (i64*)malloc((size_t)length * sizeof(i64));
    pref[0] = 1;
    for (i64 i = 0; i < length; i++)
        pref[i + 1] = (i64)((i128)pref[i] * (start + i) % mod);
    i64 inv_total = modpow(pref[length], mod - 2, mod);
    for (i64 i = length - 1; i >= 0; i--) {
        i64 x = start + i;
        invs[i] = (i64)((i128)inv_total * pref[i] % mod);
        inv_total = (i64)((i128)inv_total * x % mod);
    }
    free(pref);
    return invs;
}

long long p873_native(void) {
    i64 p = 1000000;   // 10^6
    i64 q = 10000000;  // 10^7
    i64 r = 100000000; // 10^8
    i64 mod = MOD;

    // Symmetry: ensure p <= q
    if (p > q) { i64 tmp = p; p = q; q = tmp; }

    i64 k = p + q;
    i64 tmax_transitions = 2 * p - (p == q ? 1 : 0);
    i64 tmax = tmax_transitions;
    if (r / 2 < tmax) tmax = r / 2;
    if (tmax <= 0) return 0;

    // Precompute inverses 1..p+1
    i64 *inv_small = (i64*)malloc((size_t)(p + 2) * sizeof(i64));
    inv_small[1] = 1;
    for (i64 i = 2; i <= p + 1; i++)
        inv_small[i] = (mod - (mod / i) * inv_small[mod % i] % mod) % mod;

    // Precompute C(p-1, x) for x=0..p-1
    i64 *choose_p = (i64*)malloc((size_t)p * sizeof(i64));
    choose_p[0] = 1;
    {
        i64 c = 1, n1 = p - 1;
        for (i64 x = 0; x < p - 1; x++) {
            c = (i64)((i128)c * (n1 - x) % mod);
            c = (i64)((i128)c * inv_small[x + 1] % mod);
            choose_p[x + 1] = c;
        }
    }

    // Precompute C(q-1, x) for x=0..p
    i64 *choose_q = (i64*)malloc((size_t)(p + 1) * sizeof(i64));
    choose_q[0] = 1;
    {
        i64 c = 1, n2 = q - 1;
        for (i64 x = 0; x < p; x++) {
            c = (i64)((i128)c * (n2 - x) % mod);
            c = (i64)((i128)c * inv_small[x + 1] % mod);
            choose_q[x + 1] = c;
        }
    }

    // Batch inverses for the binomial recurrence
    i64 N0 = r + k;
    i64 lowest_needed = N0 - 2 * (tmax - 1) - 1;
    i64 inv_base = lowest_needed;
    i64 inv_end = N0;
    i64 *inv_range = batch_inverses(inv_base, inv_end, mod);

    // Compute B_0 = C(r+k, k) = prod_{i=1..k} (r+i)/i
    i64 num = 1, den = 1;
    for (i64 i = 1; i <= k; i++) {
        num = (i64)((i128)num * (r + i) % mod);
        den = (i64)((i128)den * i % mod);
    }
    i64 B = (i64)((i128)num * modpow(den, mod - 2, mod) % mod);

    i64 ans = 0;
    i64 n_total = N0;
    i64 n_gap = r;

    for (i64 t = 1; t <= tmax; t++) {
        // Update B from t-1 to t
        i64 ratio = (i64)((i128)(n_gap % mod) * ((n_gap - 1) % mod) % mod);
        ratio = (i64)((i128)ratio * inv_range[n_total - inv_base] % mod);
        ratio = (i64)((i128)ratio * inv_range[n_total - 1 - inv_base] % mod);
        B = (i64)((i128)B * ratio % mod);

        n_total -= 2;
        n_gap -= 2;

        // Count AB-strings with exactly t transitions
        i64 a_runs = (t + 2) / 2;
        i64 b_runs = (t + 1) / 2;
        i64 count_t = 0;
        if (a_runs <= p && b_runs <= q)
            count_t = (i64)((i128)choose_p[a_runs - 1] * choose_q[b_runs - 1] % mod);

        i64 a_runs2 = (t + 1) / 2;
        i64 b_runs2 = (t + 2) / 2;
        if (a_runs2 <= p && b_runs2 <= q)
            count_t = (count_t + (i64)((i128)choose_p[a_runs2 - 1] * choose_q[b_runs2 - 1] % mod)) % mod;

        ans = (ans + (i64)((i128)count_t * B % mod)) % mod;
    }

    free(inv_small);
    free(choose_p);
    free(choose_q);
    free(inv_range);

    return ans;
}
