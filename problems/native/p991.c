// Project Euler 991: Equation a/(b+c) + b/(c+a) + c/(a+c) = 4
// Parameterized number theory search.
#include <stdint.h>
#include <stdio.h>
#include <math.h>

static int gcd(int a, int b) {
    while (b) { int t = a % b; a = b; b = t; }
    return a;
}

static int isqrt_ll(long long n) {
    if (n < 0) return 0;
    long long r = (long long)sqrt((double)n);
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return (int)r;
}

long long p991_native(void) {
    long long limit = 10000000;
    long long total = 0;

    // Plus branch: s = 6m^2 - n^2 + mn, n < 2m so s > 4m^2, m <= sqrt(limit/4)
    int m_max = isqrt_ll(limit / 4) + 2;
    for (int m = 1; m <= m_max; m++) {
        int n_min = isqrt_ll(3LL * m * m) + 1;
        int n_max = 2 * m - 1;
        for (int n = n_min; n <= n_max; n++) {
            if (gcd(m, n) != 1) continue;
            long long a = 4LL * m * m - (long long)n * n;
            long long c = (long long)n * n - 3LL * m * m;
            long long b = 5LL * m * m - (long long)n * n + (long long)m * n;
            long long s = a + b + c;
            if (a <= 0 || b <= 0 || c <= 0) continue;
            if (s <= limit) {
                long long count = limit / s;
                total += s * count * (count + 1) / 2;
            }
        }
    }

    // Minus branch: k = 2m - n, s = k(5m - k)
    double alpha = 2.0 + sqrt(3.0);
    double beta = (5.0 + sqrt(21.0)) / 2.0;

    for (int k = 1; ; k++) {
        int low = (int)(alpha * k) + 1;
        while ((long long)(2 * low - k) * (2 * low - k) <= 3LL * low * low)
            low++;

        int high_pos = (int)(beta * k);
        while (high_pos > 0 && !(-(long long)high_pos * high_pos + 5LL * high_pos * k - (long long)k * k > 0))
            high_pos--;

        long long high_sum = (limit + (long long)k * k) / (5 * k);
        int high = high_pos < high_sum ? high_pos : (int)high_sum;

        if (low > high_sum) break;

        for (int m = low; m <= high; m++) {
            if (gcd(m, k) != 1) continue;
            int n = 2 * m - k;
            long long a = 4LL * m * m - (long long)n * n;
            long long c = (long long)n * n - 3LL * m * m;
            long long b = 5LL * m * m - (long long)n * n - (long long)m * n;
            long long s = a + b + c;
            if (a <= 0 || b <= 0 || c <= 0) continue;
            if (s <= limit) {
                long long count = limit / s;
                total += s * count * (count + 1) / 2;
            }
        }
    }

    return total;
}
