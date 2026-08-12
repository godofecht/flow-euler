// Project Euler 993: Banana game BB(n) with eventual periodicity.
// Port of Python reference solver to C.
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define PERIOD_START 514
#define PERIOD 71
#define OFFSET 60000
#define ARR_SIZE 120001

static const int DELTA_PATTERN[PERIOD] = {
    17, -2, -8, -2, -2, -14, -2, -2, -17, -8, -5, -8, -5, -2, -2, -5, -8,
    50, -8, 23, -13, -2, 67, -5, -2, -2, -5, -8, -5, 21, 29, -11, -2, -2,
    6, -11, 31, -2, -11, 17, -2, -8, -2, -2, -14, -2, -2, -17, -8, -5, -8,
    -8, 8, -13, -5, -2, -2, -5, -2, -11, -8, -8, -5, -2, -11, -8, -8, -5,
    -2, -11, 216
};

static unsigned char bananas[ARR_SIZE];

static inline int has(int pos) {
    return bananas[pos + OFFSET];
}

static inline void set_b(int pos, int val) {
    bananas[pos + OFFSET] = (unsigned char)val;
}

static int step_state(int *pos, int *carry) {
    int p = *pos;
    int hx = has(p);
    int hx1 = has(p + 1);

    if (hx && hx1) {
        set_b(p + 1, 0);
        *pos = p - 1;
        *carry = *carry + 1;
        return 1;
    }
    if (hx && !hx1) {
        set_b(p, 0);
        *pos = p + 2;
        *carry = *carry + 1;
        return 1;
    }
    if (!hx && hx1) {
        set_b(p + 1, 0);
        set_b(p, 1);
        *pos = p + 2;
        return 1;
    }
    if (*carry >= 3) {
        set_b(p - 1, 1);
        set_b(p, 1);
        set_b(p + 1, 1);
        *pos = p - 2;
        *carry = *carry - 3;
        return 1;
    }
    return 0;
}

long long p993_native(void) {
    int limit = PERIOD_START + PERIOD;
    long long bb[586];
    bb[0] = 0;
    int pos = 0, carry = 0;
    memset(bananas, 0, sizeof(bananas));

    for (int n = 1; n <= limit; n++) {
        carry += 1;
        while (step_state(&pos, &carry)) {}
        bb[n] = pos;
    }

    long long n = 1000000000000000000LL;
    if (n <= PERIOD_START) {
        return bb[n];
    }

    long long remaining = n - PERIOD_START;
    long long whole_periods = remaining / PERIOD;
    long long tail = remaining % PERIOD;

    long long pattern_sum = 0;
    for (int i = 0; i < PERIOD; i++) pattern_sum += DELTA_PATTERN[i];

    long long tail_sum = 0;
    for (long long i = 0; i < tail; i++) tail_sum += DELTA_PATTERN[i];

    return bb[PERIOD_START] + whole_periods * pattern_sum + tail_sum;
}
