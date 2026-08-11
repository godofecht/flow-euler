#include <math.h>

/* Project Euler 724: Drone Delivery.
 *
 * E(n) = (n/2) * (H_n^2 + H_n^(2))
 *
 * For n = 10^8 we use asymptotic expansions:
 *   H_n     = log(n) + gamma + 1/(2n) - 1/(12n^2) + 1/(120n^4) - 1/(252n^6)
 *   H_n^(2) = pi^2/6 - 1/n + 1/(2n^2) - 1/(6n^3) + 1/(30n^5) - 1/(42n^7)
 *
 * Answer is round(E(10^8)) = 18128250110.
 */

#define EULER_GAMMA 0.5772156649015328606065120900824024310421

static double harmonic_asymptotic(double n) {
    double inv = 1.0 / n;
    double inv2 = inv * inv;
    double inv4 = inv2 * inv2;
    double inv6 = inv4 * inv2;
    return log(n)
        + EULER_GAMMA
        + 0.5 * inv
        - (1.0 / 12.0) * inv2
        + (1.0 / 120.0) * inv4
        - (1.0 / 252.0) * inv6;
}

static double harmonic2_asymptotic(double n) {
    double inv = 1.0 / n;
    double inv2 = inv * inv;
    double inv3 = inv2 * inv;
    double inv5 = inv3 * inv2;
    double inv7 = inv5 * inv2;
    return (M_PI * M_PI) / 6.0
        - inv
        + 0.5 * inv2
        - (1.0 / 6.0) * inv3
        + (1.0 / 30.0) * inv5
        - (1.0 / 42.0) * inv7;
}

long long p724_native(void) {
    double n = 1e8;
    double H = harmonic_asymptotic(n);
    double H2 = harmonic2_asymptotic(n);
    double e = 0.5 * n * (H * H + H2);
    return (long long)floor(e + 0.5);
}
