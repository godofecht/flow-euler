// Project Euler 853: Pisano Periods
// Sum of all n <= 10^9 with pi(n) = 120.
#include <stdint.h>
typedef long long i64;

static i64 gcd(i64 a, i64 b) {
    while (b) { i64 t = a % b; a = b; b = t; }
    return a;
}

static i64 lcm(i64 a, i64 b) {
    return a / gcd(a, b) * b;
}

static i64 primes[]     = {2,    3,   5,    11,   31,   41,   61,    2521};
static i64 pivalues[]   = {3,    8,   20,   10,   30,   40,   60,    120};
static int  mults[]     = {4,    2,   1,    1,    1,    1,    1,     1};
static int  nprimes = 8;
static i64 target = 120;
static i64 limit = 1000000000LL;

static i64 total_sum;

// Generate all products of primes[idx..nprimes-1]^k (k=0..mult) within limit.
// p = product so far, vpi = lcm of pisano periods so far.
static void gen(int idx, i64 p, i64 vpi) {
    if (idx == nprimes) {
        if (p > 1 && vpi == target)
            total_sum += p;
        return;
    }
    i64 pr = primes[idx];
    i64 pv = pivalues[idx];
    int mult = mults[idx];
    i64 cur_p = 1, cur_vpi = 1;
    for (int i = 0; i <= mult; i++) {
        if (p * cur_p < limit) {
            i64 new_vpi = (cur_p == 1) ? vpi : lcm(vpi, cur_vpi);
            gen(idx + 1, p * cur_p, new_vpi);
        }
        if (i == 0) {
            cur_p *= pr;
            cur_vpi *= pv;
        } else {
            cur_p *= pr;
            cur_vpi *= pr;
        }
    }
}

long long p853_native(void) {
    total_sum = 0;
    gen(0, 1, 1);
    return total_sum;
}
