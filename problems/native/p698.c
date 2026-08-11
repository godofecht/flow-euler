// Project Euler 698: 123 Numbers
// Find F(111111111111222333) mod 123123123.
// 123-numbers: digits only 1,2,3, and each digit count is itself a 123-number.
// Enumerate by length, count valid strings, unrank within the target length.
// Needs bignum factorials (L=38), so uses GMP.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#define MOD 123123123LL
#define TARGET_N "111111111111222333"

// 123-numbers up to 100: 1,2,3,11,12,13,21,22,23,31,32,33
// Plus 0 (meaning digit absent)
static int allowed[] = {0, 1, 2, 3, 11, 12, 13, 21, 22, 23, 31, 32, 33};
static int n_allowed = 13;
static int is_allowed[100];

static mpz_t fact[100];
static mpz_t tmp1, tmp2, tmp3;

static int is_123_number(int x) {
    if (x == 1) return 1;
    if (x <= 0) return 0;
    // Check digits
    int c1 = 0, c2 = 0, c3 = 0;
    int t = x;
    while (t > 0) {
        int d = t % 10;
        if (d != 1 && d != 2 && d != 3) return 0;
        if (d == 1) c1++;
        else if (d == 2) c2++;
        else c3++;
        t /= 10;
    }
    if (c1 && !is_123_number(c1)) return 0;
    if (c2 && !is_123_number(c2)) return 0;
    if (c3 && !is_123_number(c3)) return 0;
    return 1;
}

// Count valid strings of length L with given prefix usage (u1,u2,u3)
// A string is valid if total counts (a,b,c) with a+b+c=L, each in allowed set
static void count_completions(int L, int u1, int u2, int u3, mpz_t result) {
    int used = u1 + u2 + u3;
    int r = L - used;
    mpz_set_ui(result, 0);
    for (int ai = 0; ai < n_allowed; ai++) {
        int a = allowed[ai];
        if (a < u1) continue;
        int ra = a - u1;
        if (ra > r) continue;
        for (int bi = 0; bi < n_allowed; bi++) {
            int b = allowed[bi];
            if (b < u2) continue;
            int rb = b - u2;
            if (ra + rb > r) continue;
            int rc = r - ra - rb;
            // c = u3 + rc must be in allowed
            int c = u3 + rc;
            if (c >= 100 || !is_allowed[c]) continue;
            if (a == 0 && b == 0 && c == 0) continue;
            // multinomial: fact[r] / (fact[ra] * fact[rb] * fact[rc])
            mpz_tdiv_q(tmp1, fact[r], fact[ra]);
            mpz_tdiv_q(tmp2, tmp1, fact[rb]);
            mpz_tdiv_q(tmp3, tmp2, fact[rc]);
            mpz_add(result, result, tmp3);
        }
    }
}

long long p698_native(void) {
    // Init is_allowed
    memset(is_allowed, 0, sizeof(is_allowed));
    for (int i = 0; i < n_allowed; i++) {
        is_allowed[allowed[i]] = 1;
    }

    // Verify 123-numbers
    for (int i = 1; i < 100; i++) {
        if (is_123_number(i) && !is_allowed[i]) {
            // Found a 123-number not in our list
            fprintf(stderr, "Warning: %d is a 123-number not in allowed list\n", i);
        }
    }

    // Init factorials
    for (int i = 0; i < 100; i++) mpz_init(fact[i]);
    mpz_init_set_ui(fact[0], 1);
    for (int i = 1; i < 100; i++)
        mpz_mul_ui(fact[i], fact[i-1], i);

    mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3);

    mpz_t target, cumulative, cnt_L;
    mpz_init_set_str(target, TARGET_N, 10);
    mpz_init_set_ui(cumulative, 0);
    mpz_init(cnt_L);

    int found_L = 0;
    mpz_t rank_in_length;
    mpz_init(rank_in_length);

    for (int L = 1; L < 5000; L++) {
        // Count valid strings of length L
        count_completions(L, 0, 0, 0, cnt_L);

        mpz_add(tmp1, cumulative, cnt_L);
        if (mpz_cmp(tmp1, target) >= 0) {
            // Target is in this length
            mpz_sub(rank_in_length, target, cumulative);
            found_L = L;
            break;
        }
        mpz_set(cumulative, tmp1);
    }

    if (!found_L) {
        fprintf(stderr, "Could not find length\n");
        return -1;
    }

    // Now unrank within length found_L
    // For each position, try digit 1, 2, 3
    int u1 = 0, u2 = 0, u3 = 0;
    long long answer = 0;

    for (int pos = 0; pos < found_L; pos++) {
        for (int digit = 1; digit <= 3; digit++) {
            int nu1 = u1, nu2 = u2, nu3 = u3;
            if (digit == 1) nu1++;
            else if (digit == 2) nu2++;
            else nu3++;

            count_completions(found_L, nu1, nu2, nu3, cnt_L);

            if (mpz_cmp(rank_in_length, cnt_L) > 0) {
                mpz_sub(rank_in_length, rank_in_length, cnt_L);
            } else {
                // This digit is the one
                u1 = nu1; u2 = nu2; u3 = nu3;
                answer = (answer * 10 + digit) % MOD;
                break;
            }
        }
    }

    // Cleanup
    for (int i = 0; i < 100; i++) mpz_clear(fact[i]);
    mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3);
    mpz_clear(target); mpz_clear(cumulative); mpz_clear(cnt_L);
    mpz_clear(rank_in_length);

    return answer;
}
