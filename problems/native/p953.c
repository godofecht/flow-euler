// Project Euler 953
// Port of the Python reference solver.
// DFS over products of distinct odd primes with XOR constraint.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MOD 1000000007LL

static long long INV6;

static long long mod_pow(long long base, long long exp, long long mod) {
    long long result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = result * base % mod;
        base = base * base % mod;
        exp >>= 1;
    }
    return result;
}

static long long isqrt_ll(long long n) {
    if (n <= 0) return 0;
    long long x = (long long)sqrt((double)n);
    while (x * x > n) x--;
    while ((x + 1) * (x + 1) <= n) x++;
    return x;
}

static long long sum_sq(long long n) {
    n %= MOD;
    return n * (n + 1) % MOD * (2 * n + 1) % MOD * INV6 % MOD;
}

static long long contribution(long long kernel, long long bound, long long multiplier) {
    long long t = isqrt_ll(bound / kernel);
    return multiplier * (kernel % MOD) % MOD * sum_sq(t) % MOD;
}

static char *g_flags;
static int *g_odd_primes;
static int g_odd_prime_count;
static int g_limit;

typedef struct {
    long long prod;
    int next_idx;
    int xor_val;
    int last;
    int odd_count;
} State;

static long long branch_sum(long long bound, int target, long long multiplier) {
    if (bound <= 0) return 0;

    long long total = 0;
    if (target == 0) {
        total = contribution(1, bound, multiplier);
    }

    int plen = g_odd_prime_count;

    int stack_cap = 1 << 21;
    State *stack = (State *)malloc(stack_cap * sizeof(State));
    int sp = 0;

    for (int idx = 0; idx < plen; idx++) {
        int p = g_odd_primes[idx];
        if (p > bound / (long long)(p + 2)) break;
        if (sp >= stack_cap) {
            stack_cap *= 2;
            stack = (State *)realloc(stack, stack_cap * sizeof(State));
        }
        stack[sp].prod = p;
        stack[sp].next_idx = idx + 1;
        stack[sp].xor_val = p;
        stack[sp].last = p;
        stack[sp].odd_count = 1;
        sp++;
    }

    while (sp > 0) {
        sp--;
        State s = stack[sp];
        int next_idx = s.next_idx;
        long long prod = s.prod;
        int xor_val = s.xor_val;
        int last = s.last;
        int odd_count = s.odd_count;

        if (odd_count) {
            int cand = target ^ xor_val;
            if (cand > last && cand <= g_limit && g_flags[cand]) {
                if (prod <= bound / cand) {
                    total = (total + contribution(prod * cand, bound, multiplier)) % MOD;
                }
            }
        }

        int next_odd_count = odd_count ^ 1;
        for (int j = next_idx; j < plen; j++) {
            int q = g_odd_primes[j];
            if (prod > bound / q) break;
            long long new_prod = prod * q;
            int first = q + 2;
            if (new_prod > bound / first) break;
            if (!next_odd_count) {
                long long temp = new_prod * first;
                if (temp > bound / (first + 2)) break;
            }
            if (sp >= stack_cap) {
                stack_cap *= 2;
                stack = (State *)realloc(stack, stack_cap * sizeof(State));
            }
            stack[sp].prod = new_prod;
            stack[sp].next_idx = j + 1;
            stack[sp].xor_val = xor_val ^ q;
            stack[sp].last = q;
            stack[sp].odd_count = next_odd_count;
            sp++;
        }
    }

    free(stack);
    return total % MOD;
}

static void prime_sieve(int limit, char **flags_out, int **odd_primes_out, int *count_out) {
    char *flags = (char *)malloc(limit + 1);
    memset(flags, 1, limit + 1);
    flags[0] = flags[1] = 0;
    for (int p = 2; (long long)p * p <= limit; p++) {
        if (flags[p]) {
            for (long long j = (long long)p * p; j <= limit; j += p) {
                flags[j] = 0;
            }
        }
    }
    int count = 0;
    for (int p = 3; p <= limit; p += 2) {
        if (flags[p]) count++;
    }
    int *odd_primes = (int *)malloc(count * sizeof(int));
    int idx = 0;
    for (int p = 3; p <= limit; p += 2) {
        if (flags[p]) odd_primes[idx++] = p;
    }
    *flags_out = flags;
    *odd_primes_out = odd_primes;
    *count_out = count;
}

long long p953_native(void) {
    INV6 = mod_pow(6, MOD - 2, MOD);

    long long n = 100000000000000LL;
    int prime_limit = (int)isqrt_ll(2 * n) + 10;

    prime_sieve(prime_limit, &g_flags, &g_odd_primes, &g_odd_prime_count);
    g_limit = prime_limit;

    long long result = (branch_sum(n, 0, 1) + branch_sum(n / 2, 2, 2)) % MOD;

    free(g_flags);
    free(g_odd_primes);

    return result;
}
