// Project Euler 969: Sum S(n) for n=1..10^18 using Lagrange interpolation.
// Port of Python reference solver to C.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;
typedef __int128 i128;

static const ll MOD = 1000000007LL;

static int primes[200];
static int num_primes;

static void sieve(int n) {
    num_primes = 0;
    char *is_composite = calloc(n + 1, 1);
    for (int i = 2; i <= n; i++) {
        if (!is_composite[i]) {
            primes[num_primes++] = i;
            for (int j = i * i; j <= n; j += i)
                is_composite[j] = 1;
        }
    }
    free(is_composite);
}

static ll modpow(ll base, ll exp, ll mod) {
    base %= mod;
    if (base < 0) base += mod;
    ll result = 1;
    while (exp > 0) {
        if (exp & 1) result = result * base % mod;
        base = base * base % mod;
        exp >>= 1;
    }
    return result;
}

static ll modinv(ll x, ll mod) {
    return modpow(x % mod, mod - 2, mod);
}

// sum_{t=1}^T t^k mod MOD using Lagrange interpolation
static ll sum_pows_lagrange(ll T, int k) {
    int d = k + 1;
    if (T <= d) {
        ll s = 0;
        for (ll i = 1; i <= T; i++) {
            s = (s + modpow(i, k, MOD)) % MOD;
        }
        return s;
    }

    // y_i = sum_{t=1}^i t^k for i=0..d
    ll ys[60];
    ys[0] = 0;
    for (int i = 1; i <= d; i++) {
        ys[i] = (ys[i-1] + modpow(i, k, MOD)) % MOD;
    }

    ll fact[60];
    fact[0] = 1;
    for (int i = 1; i <= d; i++) fact[i] = fact[i-1] * i % MOD;

    ll invfact[60];
    invfact[d] = modinv(fact[d], MOD);
    for (int i = d; i >= 1; i--) invfact[i-1] = invfact[i] * i % MOD;

    ll Tmod = T % MOD;

    // pre[i] = prod_{j=0..i-1} (T - j)
    ll pre[60];
    pre[0] = 1;
    for (int i = 1; i <= d; i++) {
        pre[i] = pre[i-1] * ((Tmod - (i-1)) % MOD) % MOD;
        if (pre[i] < 0) pre[i] += MOD;
    }

    // suf[i] = prod_{j=i..d} (T - j)
    ll suf[62];
    suf[d+1] = 1;
    for (int i = d; i >= 0; i--) {
        suf[i] = suf[i+1] * ((Tmod - i) % MOD) % MOD;
        if (suf[i] < 0) suf[i] += MOD;
    }

    ll res = 0;
    for (int i = 0; i <= d; i++) {
        ll num = pre[i] * suf[i+1] % MOD;
        ll denom = fact[i] * fact[d-i] % MOD;
        if ((d - i) & 1) {
            denom = (MOD - denom) % MOD;
        }
        ll invden = modinv(denom, MOD);
        ll li = num * invden % MOD;
        res = (res + ys[i] * li) % MOD;
    }
    return res;
}

// Compute M_k = prod_{p <= k} p^{ceil(v_p(k!)/k)}
static i128 compute_Mk(int k) {
    if (k == 0) return 1;
    i128 M = 1;
    for (int pi = 0; pi < num_primes; pi++) {
        int p = primes[pi];
        if (p > k) break;
        int v = 0;
        i128 pp = p;
        while (pp <= k) {
            v += k / (int)pp;
            pp *= p;
        }
        int e = (v + k - 1) / k;  // ceil(v/k)
        if (e > 0) {
            i128 pe = 1;
            for (int i = 0; i < e; i++) pe *= p;
            M *= pe;
        }
    }
    return M;
}

static ll compute_sum_S_upto_N(ll N) {
    ll total = 0;

    // Precompute factorials mod MOD up to 200
    int MAXF = 200;
    ll fact[201];
    fact[0] = 1;
    for (int i = 1; i <= MAXF; i++) fact[i] = fact[i-1] * i % MOD;

    int k = 0;
    while (1) {
        i128 M = compute_Mk(k);
        if (M > N) break;

        ll T;
        if (k == 0) {
            T = N;
        } else {
            if (N < k) break;
            T = (N - k) / (ll)M;
            if (T <= 0) { k++; continue; }
        }

        // C_k = (-1)^k * M^k / k! mod MOD
        ll C;
        if (k == 0) {
            C = 1;
        } else {
            C = modpow((ll)(M % MOD), k, MOD) * modinv(fact[k], MOD) % MOD;
            if (k & 1) C = (MOD - C) % MOD;
        }

        ll sum_tk = sum_pows_lagrange(T, k);
        ll contrib = C * sum_tk % MOD;
        total = (total + contrib) % MOD;

        k++;
    }

    return total;
}

long long p969_native(void) {
    sieve(1000);
    return compute_sum_S_upto_N(1000000000000000000LL);
}
