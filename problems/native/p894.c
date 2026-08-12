// Project Euler 894: Spiral of Circles.
// Newton solve for (s, theta) then compute total curvilinear triangle area.
#include <math.h>
#include <stdio.h>

static double clampd(double x, double lo, double hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

static double h_func(double s, double theta, int n) {
    double sn = pow(s, n);
    double num = 1.0 + sn * sn - 2.0 * sn * cos(n * theta);
    double den = (1.0 + sn) * (1.0 + sn);
    return num / den;
}

static void f_func(double s, double theta, double *f1, double *f2) {
    double h1 = h_func(s, theta, 1);
    *f1 = h1 - h_func(s, theta, 7);
    *f2 = h1 - h_func(s, theta, 8);
}

static double objective(double s, double theta) {
    double f1, f2;
    f_func(s, theta, &f1, &f2);
    return f1 * f1 + f2 * f2;
}

static void find_initial_guess(double *s0, double *t0) {
    double best_s = 0.9, best_t = 0.8;
    double best_val = 1e300;
    for (double s = 0.75; s <= 0.99 + 1e-12; s += 0.002) {
        for (double theta = 0.05; theta <= M_PI - 0.05 + 1e-12; theta += 0.01) {
            double val = objective(s, theta);
            if (val < best_val) {
                best_val = val;
                best_s = s;
                best_t = theta;
            }
        }
    }
    *s0 = best_s;
    *t0 = best_t;
}

static void newton_solve(double s0, double t0, double *s_out, double *t_out) {
    double s = s0, theta = t0;
    for (int iter = 0; iter < 60; iter++) {
        double f1, f2;
        f_func(s, theta, &f1, &f2);
        if (fmax(fabs(f1), fabs(f2)) < 1e-15) {
            *s_out = s;
            *t_out = theta;
            return;
        }
        double ds = 1e-8, dt = 1e-8;
        double f1_sp, f2_sp, f1_sm, f2_sm, f1_tp, f2_tp, f1_tm, f2_tm;
        f_func(s + ds, theta, &f1_sp, &f2_sp);
        f_func(s - ds, theta, &f1_sm, &f2_sm);
        f_func(s, theta + dt, &f1_tp, &f2_tp);
        f_func(s, theta - dt, &f1_tm, &f2_tm);

        double a = (f1_sp - f1_sm) / (2.0 * ds);
        double b = (f1_tp - f1_tm) / (2.0 * dt);
        double c = (f2_sp - f2_sm) / (2.0 * ds);
        double d = (f2_tp - f2_tm) / (2.0 * dt);

        double det = a * d - b * c;
        if (det == 0.0 || !isfinite(det)) {
            s *= 0.999;
            theta *= 0.999;
            continue;
        }

        double delta_s = (d * f1 - b * f2) / det;
        double delta_t = (-c * f1 + a * f2) / det;

        double cur_obj = f1 * f1 + f2 * f2;
        double step = 1.0;
        int improved = 0;
        for (int ls = 0; ls < 40; ls++) {
            double ns = s - step * delta_s;
            double nt = theta - step * delta_t;
            if (!(0.0 < ns && ns < 1.0 && 0.0 < nt && nt < M_PI)) {
                step *= 0.5;
                continue;
            }
            double nobj = objective(ns, nt);
            if (nobj < cur_obj) {
                s = ns;
                theta = nt;
                improved = 1;
                break;
            }
            step *= 0.5;
        }
        if (!improved) break;
    }
    *s_out = s;
    *t_out = theta;
}

static double curvilinear_triangle_area(double r1, double r2, double r3) {
    double a = r2 + r3;
    double b = r1 + r3;
    double c = r1 + r2;
    double p = 0.5 * (a + b + c);
    double tri_sq = p * (p - a) * (p - b) * (p - c);
    if (tri_sq < 0.0) tri_sq = 0.0;
    double tri_area = sqrt(tri_sq);

    double cos1 = (b * b + c * c - a * a) / (2.0 * b * c);
    double cos2 = (a * a + c * c - b * b) / (2.0 * a * c);
    double cos3 = (a * a + b * b - c * c) / (2.0 * a * b);
    double ang1 = acos(clampd(cos1, -1.0, 1.0));
    double ang2 = acos(clampd(cos2, -1.0, 1.0));
    double ang3 = acos(clampd(cos3, -1.0, 1.0));

    double sectors = 0.5 * (r1 * r1 * ang1 + r2 * r2 * ang2 + r3 * r3 * ang3);
    return tri_area - sectors;
}

double p894_native(void) {
    double s0, t0;
    find_initial_guess(&s0, &t0);
    double s, theta;
    newton_solve(s0, t0, &s, &theta);

    double s7 = pow(s, 7);
    double s8 = s7 * s;
    double a0 = curvilinear_triangle_area(1.0, s, s8);
    double b0 = curvilinear_triangle_area(1.0, s7, s8);
    double total = (a0 + b0) / (1.0 - s * s);
    return total;
}
