// Project Euler 985: Telescoping triangles.
// Find minimum perimeter integer-sided triangle with 20 telescoping steps.
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#define EPS 1e-12

static void triangle_angles(double a, double b, double c, double *A, double *B, double *C) {
    double af = a, bf = b, cf = c;
    double cosA = (bf * bf + cf * cf - af * af) / (2.0 * bf * cf);
    double cosB = (af * af + cf * cf - bf * bf) / (2.0 * af * cf);
    double cosC = (af * af + bf * bf - cf * cf) / (2.0 * af * bf);
    if (cosA < -1.0) cosA = -1.0; if (cosA > 1.0) cosA = 1.0;
    if (cosB < -1.0) cosB = -1.0; if (cosB > 1.0) cosB = 1.0;
    if (cosC < -1.0) cosC = -1.0; if (cosC > 1.0) cosC = 1.0;
    *A = acos(cosA);
    *B = acos(cosB);
    *C = acos(cosC);
}

static int num_existing_steps(int a, int b, int c, int max_steps) {
    double A, B, C;
    triangle_angles((double)a, (double)b, (double)c, &A, &B, &C);
    int steps = 0;
    for (int i = 0; i < max_steps; i++) {
        double nA = M_PI - 2.0 * B;
        double nB = M_PI - 2.0 * C;
        double nC = M_PI - 2.0 * A;
        A = nA; B = nB; C = nC;
        if (A <= EPS || B <= EPS || C <= EPS) break;
        steps++;
    }
    return steps;
}

long long p985_native(void) {
    int target_steps = 20;
    long long best_perimeter = 0;
    int found = 0;

    for (int n = 2; n <= 5000000; n++) {
        int candidates[2][3] = {{n, n, n + 1}, {n, n + 1, n + 1}};
        for (int ci = 0; ci < 2; ci++) {
            int a = candidates[ci][0], b = candidates[ci][1], c = candidates[ci][2];
            int steps = num_existing_steps(a, b, c, target_steps + 2);
            if (steps == target_steps) {
                int p = a + b + c;
                if (!found || p < best_perimeter) {
                    best_perimeter = p;
                    found = 1;
                }
            }
        }
        if (found && 3 * (long long)(n + 1) + 1 > best_perimeter) {
            break;
        }
    }
    return best_perimeter;
}
