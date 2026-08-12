#include <stdint.h>
#include <math.h>

typedef int64_t i64;
typedef unsigned long long u64;

static const i64 MOD = 1000000007LL;

static u64 isqrt_u64(u64 n) {
    if (n < 2) return n;
    u64 x = (u64)sqrtl((long double)n);
    while (x * x > n) x--;
    while ((x + 1) * (x + 1) <= n) x++;
    return x;
}

/* Closed form for the "easy" prefix sum:
 * (m^3 + 15 m^2 - 52 m + 36) / 6  mod MOD   (m >= 3)
 */
static i64 easy_prefix(i64 m) {
    if (m < 3) return 0;
    /* Work modulo 6*MOD to keep division exact, then divide by 6. */
    i64 mm = MOD * 6;
    i64 m3 = ((__int128)m * m % mm) * m % mm;
    i64 m2 = (__int128)m * m % mm;
    i64 num = (m3 + 15 * m2 % mm + (mm - 52 * m % mm) + 36) % mm;
    return num / 6;
}

/* For fixed (x,y), compute n_min and whether D is a perfect square.
 * n_min = 2xy + x + y + 1 + ceil(2*sqrt(D)), D = x*y*(x+1)*(y+1).
 * ceil(2*sqrt(D)) = ceil(sqrt(4D)).  If 4D is a perfect square, sq=1.
 */
static i64 min_n_and_square(i64 x, i64 y, int *sq) {
    u64 four_d = (u64)4 * (u64)(x * (x + 1)) * (u64)(y * (y + 1));
    u64 r = isqrt_u64(four_d);
    if (r * r == four_d) {
        *sq = 1;
        return 2 * x * y + x + y + 1 + (i64)r;
    }
    *sq = 0;
    return 2 * x * y + x + y + 1 + (i64)(r + 1);
}

/* Binary search for max y >= x with n_min(x,y) <= m. */
static i64 y_max_for_x(i64 m, i64 x) {
    if (4 * x > m - 1) return x - 1;
    i64 hi = (m - 1) / (4 * x) + 2;
    if (hi < x) hi = x;
    i64 lo = x, ok = x - 1;
    while (lo <= hi) {
        i64 mid = (lo + hi) / 2;
        int sq;
        i64 n_min = min_n_and_square(x, mid, &sq);
        if (n_min <= m) { ok = mid; lo = mid + 1; }
        else hi = mid - 1;
    }
    return ok;
}

static i64 hard_prefix(i64 m) {
    if (m < 3) return 0;
    i64 k = (m - 1) / 4;
    if (k <= 0) return 0;
    i64 x_max = (i64)isqrt_u64((u64)k);

    i64 total = 0;
    i64 cutoff = (i64)MOD << 20;

    for (i64 x = 1; x <= x_max; x++) {
        i64 ymax = y_max_for_x(m, x);
        if (ymax < x) continue;

        i64 A = x * (x + 1);
        i64 y = x;
        i64 yy1 = y * (y + 1);
        i64 two_xy = 2 * x * y;

        while (y <= ymax) {
            u64 four_d = (u64)(A * yy1) << 2;
            u64 r = isqrt_u64(four_d);
            int sq = (r * r == four_d) ? 1 : 0;
            i64 ceil2 = sq ? (i64)r : (i64)(r + 1);
            i64 n_min = two_xy + x + y + 1 + ceil2;

            if (n_min <= m) {
                i64 cnt = 2 * (m - n_min + 1) - sq;
                i64 add = (x == y) ? cnt : cnt << 1;
                total += add;
                if (total >= cutoff) total %= MOD;
            }

            yy1 += (y << 1) + 2;
            y++;
            two_xy += x << 1;
        }
    }
    return total % MOD;
}

long long p747_native(void) {
    i64 m = 100000000LL;
    i64 easy = easy_prefix(m);
    i64 hard = hard_prefix(m);
    i64 result = (easy + 3 * hard) % MOD;
    return (long long)result;
}
