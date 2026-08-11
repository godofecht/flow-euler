/* Project Euler 855 - Delphi Paper.
 *
 * Closed form:
 *   S(a,b) = (a!)^b * (b!)^a / ((ab)!)^2
 *
 * Computed in log-space to avoid overflow, then exponentiated.
 * With a=5, b=8 the value is ~6.88e-57, well within double range.
 */

#include <math.h>

double p855_native(void) {
    int a = 5;
    int b = 8;
    double log_s = (double)b * lgamma((double)(a + 1))
                 + (double)a * lgamma((double)(b + 1))
                 - 2.0 * lgamma((double)(a * b + 1));
    return exp(log_s);
}
