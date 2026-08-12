// Project Euler 968: Counting via determinants, adjugates, and modular arithmetic.
// Port of Python reference solver to C.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;

static const ll MOD = 1000000007LL;

static const int primes[5] = {2, 3, 5, 7, 11};

static const int h[15][5] = {
    {1, 1, 0, 0, 0},
    {1, 0, 1, 0, 0},
    {1, 0, 0, 1, 0},
    {1, 0, 0, 0, 1},
    {0, 1, 1, 0, 0},
    {0, 1, 0, 1, 0},
    {0, 1, 0, 0, 1},
    {0, 0, 1, 1, 0},
    {0, 0, 1, 0, 1},
    {0, 0, 0, 1, 1},
    {-1, 0, 0, 0, 0},
    {0, -1, 0, 0, 0},
    {0, 0, -1, 0, 0},
    {0, 0, 0, -1, 0},
    {0, 0, 0, 0, -1},
};

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

// Determinant of n x n matrix (n <= 5), integer arithmetic
static ll det_matrix(int n, ll M[5][5]) {
    if (n == 1) return M[0][0];
    if (n == 2) return M[0][0]*M[1][1] - M[0][1]*M[1][0];
    ll s = 0;
    for (int j = 0; j < n; j++) {
        if (M[0][j] == 0) continue;
        ll sub[5][5];
        for (int i = 1; i < n; i++) {
            int col = 0;
            for (int k = 0; k < n; k++) {
                if (k == j) continue;
                sub[i-1][col] = M[i][k];
                col++;
            }
        }
        ll cof = det_matrix(n-1, sub);
        if (j & 1) s -= M[0][j] * cof;
        else s += M[0][j] * cof;
    }
    return s;
}

// Adjugate of n x n matrix
static void adj_matrix(int n, ll M[5][5], ll adj[5][5]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            ll sub[5][5];
            int row = 0;
            for (int k = 0; k < n; k++) {
                if (k == j) continue;
                int col = 0;
                for (int l = 0; l < n; l++) {
                    if (l == i) continue;
                    sub[row][col] = M[k][l];
                    col++;
                }
                row++;
            }
            ll cof = det_matrix(n-1, sub);
            if ((i+j) & 1) adj[i][j] = -cof;
            else adj[i][j] = cof;
        }
    }
}

// w(v): product of primes^v[j] mod m, handling negative exponents
static ll w_func(ll v[5]) {
    ll z = 1;
    for (int j = 0; j < 5; j++) {
        ll e = v[j];
        if (e >= 0) {
            z = z * modpow(primes[j], e, MOD) % MOD;
        } else {
            ll inv_base = modinv(modpow(primes[j], -e, MOD), MOD);
            z = z * inv_base % MOD;
        }
    }
    return z;
}

// f(v): product of primes^v[j] mod m, return 0 if any exponent < 0
static ll f_func(ll v[5]) {
    for (int j = 0; j < 5; j++) {
        if (v[j] < 0) return 0;
    }
    ll z = 1;
    for (int j = 0; j < 5; j++) {
        z = z * modpow(primes[j], v[j], MOD) % MOD;
    }
    return z;
}

// Stored combination data
typedef struct {
    int I[5];
    ll det;
    ll adj[5][5];
    ll r[5][5];  // signed columns of adj
    ll den;      // modular denominator product
} Combo;

static Combo combos[3003];
static int num_combos;

static void precompute_combos() {
    num_combos = 0;
    // Generate all C(15,5) combinations
    for (int a = 0; a < 15; a++)
    for (int b = a+1; b < 15; b++)
    for (int c = b+1; c < 15; c++)
    for (int d = c+1; d < 15; d++)
    for (int e = d+1; e < 15; e++) {
        int idx[5] = {a, b, c, d, e};
        ll M[5][5];
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
                M[i][j] = h[idx[i]][j];

        ll det = det_matrix(5, M);
        if (det != 1 && det != -1 && det != 2 && det != -2) continue;

        ll adj[5][5];
        adj_matrix(5, M, adj);

        // Build r: columns of adj, possibly negated
        ll r[5][5];
        for (int i = 0; i < 5; i++) {
            ll v[5];
            for (int j = 0; j < 5; j++) v[j] = adj[j][i];
            ll dot = 0;
            for (int k = 0; k < 5; k++) dot += M[i][k] * v[k];
            if (dot > 0) {
                for (int j = 0; j < 5; j++) r[i][j] = -v[j];
            } else {
                for (int j = 0; j < 5; j++) r[i][j] = v[j];
            }
        }

        // den = product of (1 - w(r[i]))^(-1) mod m
        ll den = 1;
        for (int i = 0; i < 5; i++) {
            ll wi = w_func(r[i]);
            ll term = (1 - wi) % MOD;
            if (term < 0) term += MOD;
            den = den * modinv(term, MOD) % MOD;
        }

        Combo *cb = &combos[num_combos++];
        for (int i = 0; i < 5; i++) cb->I[i] = idx[i];
        cb->det = det;
        memcpy(cb->adj, adj, sizeof(adj));
        memcpy(cb->r, r, sizeof(r));
        cb->den = den;
    }
}

// P function: compute P(L) where L has 10 values
static ll P_func(ll L[10], int s_flag) {
    ll b[15];
    for (int i = 0; i < 10; i++) b[i] = L[i];
    for (int i = 10; i < 15; i++) b[i] = 0;

    ll t = 0;
    for (int ci = 0; ci < num_combos; ci++) {
        Combo *cb = &combos[ci];
        ll det = cb->det;

        ll bn[5];
        for (int i = 0; i < 5; i++) bn[i] = b[cb->I[i]];

        ll xn[5];
        for (int i = 0; i < 5; i++) {
            ll sum = 0;
            for (int j = 0; j < 5; j++) sum += cb->adj[i][j] * bn[j];
            xn[i] = sum;
        }

        int ok = 1;
        int eq = 0;
        if (det > 0) {
            for (int idx = 0; idx < 15 && ok; idx++) {
                ll v = 0;
                for (int j = 0; j < 5; j++) v += h[idx][j] * xn[j];
                ll u = b[idx] * det;
                if (v > u) { ok = 0; break; }
                if (v == u) eq++;
            }
        } else {
            for (int idx = 0; idx < 15 && ok; idx++) {
                ll v = 0;
                for (int j = 0; j < 5; j++) v += h[idx][j] * xn[j];
                ll u = b[idx] * det;
                if (v < u) { ok = 0; break; }
                if (v == u) eq++;
            }
        }

        if (!ok || (s_flag && eq != 5)) continue;

        ll da = det < 0 ? -det : det;
        if (da == 1) {
            ll exp[5];
            for (int i = 0; i < 5; i++) exp[i] = xn[i] / det;  // exact since det=±1
            t = (t + f_func(exp) * cb->den) % MOD;
        } else {
            // da == 2
            ll dr[5][5];
            for (int i = 0; i < 5; i++)
                for (int j = 0; j < 5; j++)
                    dr[i][j] = det * cb->r[i][j];

            ll tw = da * 2;  // 4
            ll num = 0;
            for (int mask = 0; mask < 32; mask++) {
                ll T[5];
                for (int j = 0; j < 5; j++) T[j] = 2 * xn[j];
                for (int i = 0; i < 5; i++) {
                    if (mask & (1 << i)) {
                        for (int j = 0; j < 5; j++) T[j] += dr[i][j];
                    }
                }
                int valid = 1;
                for (int j = 0; j < 5; j++) {
                    if (T[j] % tw != 0) { valid = 0; break; }
                }
                if (!valid) continue;

                ll exp[5];
                for (int j = 0; j < 5; j++) exp[j] = T[j] / (2 * det);

                int all_ok = 1;
                for (int row = 0; row < 15 && all_ok; row++) {
                    ll v = 0;
                    for (int k = 0; k < 5; k++) v += h[row][k] * exp[k];
                    if (v > b[row]) { all_ok = 0; break; }
                }
                if (all_ok) {
                    num = (num + f_func(exp)) % MOD;
                }
            }
            t = (t + num * cb->den) % MOD;
        }
    }
    return t;
}

static ll B_sequence(int n, ll *a) {
    a[0] = 1;
    a[1] = 7;
    for (int i = 2; i < n; i++) {
        a[i] = (7 * a[i-1] + a[i-2] * a[i-2]) % MOD;
    }
    return 0;
}

long long p968_native(void) {
    precompute_combos();

    int n = 100;
    ll A[1005];
    B_sequence(10 * n + 5, A);

    ll total = 0;
    for (int i = 0; i < n; i++) {
        ll group[10];
        for (int j = 0; j < 10; j++) group[j] = A[10 * i + j];
        total = (total + P_func(group, 1)) % MOD;
    }
    return total;
}
