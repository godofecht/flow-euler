// Project Euler 621: Sum of three triangular numbers
// G(n) = r_3(8n+3) / 8
// Uses class number formula: G(n) = 3 * h(D) * S / w_div2
// where D = -n0, N = n0 * f^2, S = sum_{d|f} mu(d)*(D/d)*sigma(f/d)
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Number theory helpers
static long long isqrt_ll(long long n) {
    if (n < 0) return 0;
    long long x = n, y = (n + 1) / 2;
    while (y < x) { x = y; y = (x + n / x) / 2; }
    return x;
}

static long long gcd_ll(long long a, long long b) {
    while (b) { long long t = a % b; a = b; b = t; }
    return a;
}

static long long mul_mod_ll(long long a, long long b, long long m) {
    return (a % m) * (b % m) % m;
}

static long long pow_mod_ll(long long a, long long e, long long m) {
    long long r = 1, b = a % m;
    while (e > 0) {
        if (e & 1) r = r * b % m;
        b = b * b % m;
        e >>= 1;
    }
    return r;
}

// Miller-Rabin (deterministic for 64-bit)
static int is_prime_ll(long long n) {
    if (n < 2) return 0;
    int small[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++) {
        if (n == small[i]) return 1;
        if (n % small[i] == 0) return 0;
    }
    long long d = n - 1;
    int s = 0;
    while ((d & 1) == 0) { d >>= 1; s++; }
    int bases[] = {2,3,5,7,11,13,17};
    for (int i = 0; i < 7; i++) {
        long long a = bases[i];
        if (a % n == 0) continue;
        long long x = pow_mod_ll(a, d, n);
        if (x == 1 || x == n - 1) continue;
        int found = 0;
        for (int j = 0; j < s - 1; j++) {
            x = x * x % n;
            if (x == n - 1) { found = 1; break; }
        }
        if (!found) return 0;
    }
    return 1;
}

// Pollard's rho
static long long pollard_rho(long long n) {
    if (n % 2 == 0) return 2;
    if (n % 3 == 0) return 3;
    if (is_prime_ll(n)) return n;
    long long c = 1, x = 2, y = 2, d = 1;
    while (d == 1) {
        x = (x * x + c) % n;
        y = (y * y + c) % n;
        y = (y * y + c) % n;
        d = gcd_ll(x > y ? x - y : y - x, n);
    }
    if (d == n) {
        // retry with different c
        for (c = 2; c < 100; c++) {
            x = y = 2; d = 1;
            while (d == 1) {
                x = (x * x + c) % n;
                y = (y * y + c) % n;
                y = (y * y + c) % n;
                d = gcd_ll(x > y ? x - y : y - x, n);
            }
            if (d != n) return d;
        }
    }
    return d;
}

// Factorize n into primes[] and exps[]
static int factorize_ll(long long n, long long *primes, int *exps) {
    int cnt = 0;
    while (n > 1) {
        if (is_prime_ll(n)) {
            int found = 0;
            for (int i = 0; i < cnt; i++) {
                if (primes[i] == n) { exps[i]++; found = 1; break; }
            }
            if (!found) {
                primes[cnt] = n;
                exps[cnt] = 1;
                cnt++;
            }
            break;
        }
        long long d = pollard_rho(n);
        // find smallest factor
        while (!is_prime_ll(d)) d = pollard_rho(d);
        // add to list
        int found = 0;
        for (int i = 0; i < cnt; i++) {
            if (primes[i] == d) { exps[i]++; found = 1; break; }
        }
        if (!found) {
            primes[cnt] = d;
            exps[cnt] = 1;
            cnt++;
        }
        n /= d;
        done:;
    }
    return cnt;
}

// Legendre symbol (a/p) for odd prime p
static int legendre(long long a, long long p) {
    a %= p;
    if (a < 0) a += p;
    if (a == 0) return 0;
    long long t = pow_mod_ll(a, (p - 1) / 2, p);
    return (t == p - 1) ? -1 : 1;
}

// Tonelli-Shanks: sqrt of n mod p (odd prime)
static long long tonelli_shanks(long long n, long long p) {
    n %= p;
    if (n < 0) n += p;
    if (n == 0) return 0;
    if (legendre(n, p) != 1) return -1;
    if (p % 4 == 3) return pow_mod_ll(n, (p + 1) / 4, p);

    long long q = p - 1;
    int s = 0;
    while ((q & 1) == 0) { q >>= 1; s++; }

    long long z = 2;
    while (legendre(z, p) != -1) z++;

    long long c = pow_mod_ll(z, q, p);
    long long r = pow_mod_ll(n, (q + 1) / 2, p);
    long long t = pow_mod_ll(n, q, p);
    int m = s;

    while (t != 1) {
        int i = 1;
        long long t2 = t * t % p;
        while (i < m && t2 != 1) {
            t2 = t2 * t2 % p;
            i++;
        }
        long long b = pow_mod_ll(c, 1LL << (m - i - 1), p);
        r = r * b % p;
        t = t * b % p * b % p;
        c = b * b % p;
        m = i;
    }
    return r;
}

// Hensel lift: lift root r mod p to root mod p^e
static long long hensel_lift(long long n, long long p, int e, long long r) {
    long long pe = p;
    long long r_mod = r % p;
    // Normalize n to be positive mod pe (final pe = p^e)
    long long pe_final = 1;
    for (int i = 0; i < e; i++) pe_final *= p;
    n %= pe_final;
    if (n < 0) n += pe_final;
    for (int k = 1; k < e; k++) {
        long long diff = r_mod * r_mod - n;
        // Floor division (Python-style) for diff / pe
        long long q = diff / pe;
        if (diff % pe != 0 && ((diff < 0) != (pe < 0))) q -= 1;
        long long rhs = q % p;
        if (rhs < 0) rhs += p;
        long long inv = pow_mod_ll((2 * r_mod) % p, p - 2, p);
        long long t = (p - rhs) * inv % p;
        r_mod = r_mod + t * pe;
        pe *= p;
        r_mod %= pe;
    }
    return r_mod;
}

// Roots of x^2 ≡ n (mod p^e) for odd prime p
// Returns number of roots, fills roots[]
static int roots_mod_pe(long long n, long long p, int e, long long *roots) {
    long long pe = 1;
    for (int i = 0; i < e; i++) pe *= p;
    n %= pe;
    if (n < 0) n += pe;

    if (n % p == 0) {
        if (n == 0) {
            long long step = 1;
            for (int i = 0; i < (e + 1) / 2; i++) step *= p;
            int cnt = 0;
            for (long long x = 0; x < pe; x += step) {
                roots[cnt++] = x;
            }
            return cnt;
        }
        if (e == 1) { roots[0] = 0; return 1; }
        return 0;
    }

    long long r = tonelli_shanks(n, p);
    if (r == -1) return 0;
    if (e > 1) r = hensel_lift(n, p, e, r);
    long long r2 = pe - r;
    if (r2 == r) { roots[0] = r; return 1; }
    roots[0] = r; roots[1] = r2;
    return 2;
}

// Extended Euclidean algorithm: returns gcd(a, b), sets x and y such that a*x + b*y = gcd
static long long ext_gcd(long long a, long long b, long long *x, long long *y) {
    if (b == 0) { *x = 1; *y = 0; return a; }
    long long x1, y1;
    long long g = ext_gcd(b, a % b, &x1, &y1);
    *x = y1;
    *y = x1 - (a / b) * y1;
    return g;
}

// Modular inverse of a mod m (m not necessarily prime, but gcd(a,m)=1)
static long long mod_inv_general(long long a, long long m) {
    long long x, y;
    a %= m;
    if (a < 0) a += m;
    ext_gcd(a, m, &x, &y);
    x %= m;
    if (x < 0) x += m;
    return x;
}

// CRT: combine x≡a1 (mod m1), x≡a2 (mod m2)
static void crt(long long a1, long long m1, long long a2, long long m2,
                long long *out_x, long long *out_m) {
    if (m1 == 1) { *out_x = a2; *out_m = m2; return; }
    if (m2 == 1) { *out_x = a1; *out_m = m1; return; }
    long long inv = mod_inv_general(m1 % m2, m2);
    long long t = ((a2 - a1) % m2 + m2) % m2 * inv % m2;
    *out_x = a1 + m1 * t;
    *out_m = m1 * m2;
}

// Sieve smallest prime factor
static int *spf_arr;
static void sieve_spf(int n) {
    spf_arr = (int*)malloc((n + 1) * sizeof(int));
    for (int i = 0; i <= n; i++) spf_arr[i] = i;
    for (int i = 2; i * i <= n; i++) {
        if (spf_arr[i] == i) {
            for (int j = i * i; j <= n; j += i) {
                if (spf_arr[j] == j) spf_arr[j] = i;
            }
        }
    }
}

// Factorize small number using spf
static int factorize_small(int x, int *primes, int *exps) {
    int cnt = 0;
    while (x > 1) {
        int p = spf_arr[x];
        int e = 0;
        while (x % p == 0) { x /= p; e++; }
        primes[cnt] = p;
        exps[cnt] = e;
        cnt++;
    }
    return cnt;
}

// Class number h(D) for D < 0, D ≡ 5 mod 8
static long long class_number(long long D) {
    long long absD = -D;
    long long amax = isqrt_ll(absD / 3);

    sieve_spf((int)amax);

    long long h = 0;

    // For D ≡ 5 mod 8, only odd a values
    for (long long a = 1; a <= amax; a += 2) {
        int fac_p[20], fac_e[20];
        int nfac = factorize_small((int)a, fac_p, fac_e);

        // Check: if p|D and p^2|a, skip
        int bad = 0;
        for (int i = 0; i < nfac; i++) {
            if (fac_e[i] >= 2 && D % fac_p[i] == 0) { bad = 1; break; }
        }
        if (bad) continue;

        // Build roots mod a via CRT
        long long roots[1024];
        int n_roots = 1;
        roots[0] = 0;
        long long mod = 1;
        int ok = 1;

        for (int i = 0; i < nfac; i++) {
            long long p = fac_p[i];
            int e = fac_e[i];
            long long pe = 1;
            for (int j = 0; j < e; j++) pe *= p;

            long long rset[8];
            int nr = roots_mod_pe(D, p, e, rset);
            if (nr == 0) { ok = 0; break; }

            long long new_roots[1024];
            int new_nr = 0;
            for (int j = 0; j < n_roots; j++) {
                for (int k = 0; k < nr; k++) {
                    long long x, m;
                    crt(roots[j], mod, rset[k], pe, &x, &m);
                    new_roots[new_nr++] = x;
                }
            }
            memcpy(roots, new_roots, new_nr * sizeof(long long));
            n_roots = new_nr;
            mod *= pe;
        }
        if (!ok) continue;

        for (int ri = 0; ri < n_roots; ri++) {
            long long r = roots[ri];
            long long b;
            if (r == 0) {
                b = a;
            } else {
                b = (r & 1) ? r : (r - (long long)a);
            }

            if (b < 0) b = -b;
            if (b > a) continue;

            long long num = b * b - D;
            long long den = 4 * a;
            if (num % den != 0) continue;
            long long c = num / den;
            if (a > c) continue;

            // reduced boundary rule
            // recompute signed b
            long long b_signed;
            if (r == 0) b_signed = a;
            else b_signed = (r & 1) ? r : (r - (long long)a);

            if ((b_signed < 0 || b_signed == a) && (b_signed == a || a == c) && b_signed < 0) continue;
            // Simplified: if (|b|==a or a==c) and b<0, skip
            long long abs_b = b_signed < 0 ? -b_signed : b_signed;
            if ((abs_b == a || a == c) && b_signed < 0) continue;

            h++;
        }
    }

    free(spf_arr);
    return h;
}

long long p621_native(void) {
    long long n = 17526LL * 1000000000LL;
    long long N = 8 * n + 3;

    // Factor N
    long long primes[20];
    int exps[20];
    int nfac = factorize_ll(N, primes, exps);

    // N = n0 * f^2, n0 squarefree
    long long n0 = 1, f = 1;
    for (int i = 0; i < nfac; i++) {
        if (exps[i] & 1) n0 *= primes[i];
        for (int j = 0; j < exps[i] / 2; j++) f *= primes[i];
    }

    long long D = -n0;
    long long h = class_number(D);

    // w_div2
    long long w_div2;
    if (D == -3) w_div2 = 3;
    else if (D == -4) w_div2 = 2;
    else w_div2 = 1;

    // S = sum_{d|f} mu(d) * (D/d) * sigma(f/d)
    // f = 81 = 3^4, primes of f = [3], exps = [4]
    long long f_primes[20];
    int f_exps[20];
    int nffac = factorize_ll(f, f_primes, f_exps);

    long long S = 0;
    int nprimes = nffac;
    for (int mask = 0; mask < (1 << nprimes); mask++) {
        int bits = 0;
        long long jac = 1;
        long long sig = 1;
        int valid = 1;
        for (int i = 0; i < nprimes; i++) {
            long long p = f_primes[i];
            int e = f_exps[i];
            if ((mask >> i) & 1) {
                bits++;
                jac *= legendre(D, p);
                e--;
            }
            if (e < 0) { sig = 0; valid = 0; break; }
            // sigma(p^e) = (p^(e+1)-1)/(p-1)
            long long pe = 1;
            for (int j = 0; j <= e; j++) pe *= p;
            sig *= (pe - 1) / (p - 1);
        }
        if (!valid || sig == 0) continue;
        long long sign = (bits & 1) ? -1 : 1;
        S += sign * jac * sig;
    }

    // G(n) = 3 * h * S / w_div2
    long long result = 3 * h * S / w_div2;
    return result;
}
