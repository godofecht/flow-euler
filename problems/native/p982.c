// Project Euler 982: Dice game optimal strategy via linear programming.
// Direct port of Python reference solver (two-phase simplex).
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EPS 1e-9

// Dynamic tableau: rows are double* arrays.
static double **tab;
static int *basis;
static int M, N;  // M rows (including objective), N cols (including RHS)

static void pivot_op(int row, int col) {
    double pv = tab[row][col];
    double inv = 1.0 / pv;
    for (int j = 0; j < N; j++) tab[row][j] *= inv;
    for (int i = 0; i < M; i++) {
        if (i == row) continue;
        double f = tab[i][col];
        if (fabs(f) > EPS) {
            for (int j = 0; j < N; j++) tab[i][j] -= f * tab[row][j];
        }
    }
    basis[row] = col;
}

static int simplex_max() {
    int m = M - 1, n = N - 1;
    for (int iter = 0; iter < 2000000; iter++) {
        int entering = -1;
        for (int j = 0; j < n; j++) {
            if (tab[m][j] < -EPS) { entering = j; break; }
        }
        if (entering == -1) return 1;
        double min_ratio = 1e18;
        int leaving = -1;
        for (int i = 0; i < m; i++) {
            double a = tab[i][entering];
            if (a > EPS) {
                double ratio = tab[i][n] / a;
                if (ratio < min_ratio - EPS) { min_ratio = ratio; leaving = i; }
            }
        }
        if (leaving == -1) return 0;
        pivot_op(leaving, entering);
    }
    return -1;
}

static void set_objective(double *c) {
    int m = M - 1, n = N - 1;
    for (int j = 0; j < N; j++) tab[m][j] = 0.0;
    for (int j = 0; j < n; j++) tab[m][j] = -c[j];
    for (int i = 0; i < m; i++) {
        double cb = c[basis[i]];
        if (fabs(cb) > EPS) {
            for (int j = 0; j < N; j++) tab[m][j] += cb * tab[i][j];
        }
    }
}

typedef struct { double *coeffs; char sense; double rhs; } Constraint;

static int *art_flags;
static int num_art;
static int *art_indices;

static void build_tableau(int n_vars, Constraint *cons, int nc) {
    // First pass: count total variables.
    int n_total = n_vars;
    for (int ci = 0; ci < nc; ci++) {
        char s = cons[ci].sense;
        if (s == '<') n_total += 1;
        else if (s == '>') n_total += 2;
        else n_total += 1;
    }
    M = nc + 1;  // +1 for objective row (added later)
    N = n_total + 1;  // +1 for RHS
    tab = malloc(M * sizeof(double *));
    for (int i = 0; i < M; i++) tab[i] = calloc(N, sizeof(double));
    basis = malloc(M * sizeof(int));
    art_flags = calloc(n_total, sizeof(int));
    art_indices = malloc(n_total * sizeof(int));
    num_art = 0;

    int cur = n_vars;
    for (int ci = 0; ci < nc; ci++) {
        double b = cons[ci].rhs;
        char sense = cons[ci].sense;
        if (b < 0) {
            for (int j = 0; j < n_vars; j++) tab[ci][j] = -cons[ci].coeffs[j];
            b = -b;
            if (sense == '<') sense = '>';
            else if (sense == '>') sense = '<';
        } else {
            for (int j = 0; j < n_vars; j++) tab[ci][j] = cons[ci].coeffs[j];
        }
        if (sense == '<') {
            tab[ci][cur] = 1.0; basis[ci] = cur; cur++;
        } else if (sense == '>') {
            tab[ci][cur] = -1.0; cur++;
            tab[ci][cur] = 1.0; basis[ci] = cur;
            art_flags[cur] = 1; art_indices[num_art++] = cur; cur++;
        } else {  // '='
            tab[ci][cur] = 1.0; basis[ci] = cur;
            art_flags[cur] = 1; art_indices[num_art++] = cur; cur++;
        }
        tab[ci][N - 1] = b;  // RHS at last column
    }
    M = nc;  // currently only constraint rows; objective added by caller
}

static void remove_artificial() {
    int m = M;  // constraint rows
    int n = N - 1;  // variables

    // Try to pivot artificial variables out of basis.
    int i = 0;
    while (i < m) {
        if (art_flags[basis[i]]) {
            int pc = -1;
            for (int j = 0; j < n; j++) {
                if (art_flags[j]) continue;
                if (fabs(tab[i][j]) > EPS) { pc = j; break; }
            }
            if (pc != -1) {
                pivot_op(i, pc);
                i++;
            } else {
                // Redundant row; remove it.
                free(tab[i]);
                for (int r = i; r < m - 1; r++) {
                    tab[r] = tab[r + 1];
                    basis[r] = basis[r + 1];
                }
                tab[m - 1] = NULL;
                m--;
            }
        } else {
            i++;
        }
    }

    // Remove artificial columns.
    int n_keep = 0;
    int *keep = malloc(n * sizeof(int));
    for (int j = 0; j < n; j++) {
        if (!art_flags[j]) keep[n_keep++] = j;
    }

    // Build mapping from old to new column index.
    int *mapping = malloc(n * sizeof(int));
    int new_idx = 0;
    for (int j = 0; j < n; j++) {
        if (!art_flags[j]) mapping[j] = new_idx++;
        else mapping[j] = -1;
    }

    // Rebuild tableau with only kept columns + RHS.
    int new_N = n_keep + 1;
    for (int r = 0; r < m; r++) {
        double *new_row = calloc(new_N, sizeof(double));
        for (int j = 0; j < n_keep; j++)
            new_row[j] = tab[r][keep[j]];
        new_row[n_keep] = tab[r][n];  // RHS
        free(tab[r]);
        tab[r] = new_row;
        basis[r] = mapping[basis[r]];
    }

    N = new_N;
    M = m;
    free(keep);
    free(mapping);
}

static double solve_lp(int n_vars, Constraint *cons, int nc, double *obj) {
    build_tableau(n_vars, cons, nc);

    if (num_art > 0) {
        // Phase I: maximize -sum(artificial).
        int nt = N - 1;
        double *c1 = calloc(nt, sizeof(double));
        for (int j = 0; j < nt; j++) if (art_flags[j]) c1[j] = -1.0;

        // Add objective row.
        int obj_row = M;
        M++;
        tab = realloc(tab, M * sizeof(double *));
        tab[obj_row] = calloc(N, sizeof(double));
        set_objective(c1);
        simplex_max();

        if (tab[obj_row][N - 1] < -1e-7) {
            // Infeasible.
            for (int i = 0; i < M; i++) free(tab[i]);
            free(tab); free(basis); free(art_flags); free(art_indices); free(c1);
            return -1e18;
        }

        // Remove objective row.
        free(tab[obj_row]);
        M--;

        // Remove artificial variables.
        remove_artificial();
        free(c1);
    }

    // Phase II: minimize sum(obj[j] * x[j]) = maximize -sum(obj[j] * x[j]).
    int nt2 = N - 1;
    double *c2 = calloc(nt2, sizeof(double));
    for (int j = 0; j < n_vars; j++) c2[j] = -obj[j];

    // Add objective row.
    int obj_row = M;
    M++;
    tab = realloc(tab, M * sizeof(double *));
    tab[obj_row] = calloc(N, sizeof(double));
    set_objective(c2);
    simplex_max();

    double result = -tab[obj_row][N - 1];

    for (int i = 0; i < M; i++) free(tab[i]);
    free(tab); free(basis); free(art_flags); free(art_indices); free(c2);
    return result;
}

static double build_and_solve(int num_dice) {
    int values[] = {1, 2, 3, 4, 5, 6};
    int num_states = 1;
    for (int i = 0; i < num_dice; i++) num_states *= 6;

    int (*states)[3] = malloc(num_states * sizeof(*states));
    int s_idx = 0;
    for (int a = 0; a < 6; a++)
    for (int b = 0; b < 6; b++)
    for (int c = 0; c < 6; c++) {
        states[s_idx][0] = values[a];
        states[s_idx][1] = values[b];
        states[s_idx][2] = values[c];
        s_idx++;
    }

    int sig_arr[100][2];
    int num_signals = 0;
    for (int a = 0; a < 6; a++) {
        for (int b = a; b < 6; b++) {
            sig_arr[num_signals][0] = values[a];
            sig_arr[num_signals][1] = values[b];
            num_signals++;
        }
    }

    int num_x = num_states * num_dice;
    int num_z = num_signals;
    int n_vars = num_x + num_z;

    int n_constraints = num_states + 2 * num_signals;
    Constraint *constraints = malloc(n_constraints * sizeof(Constraint));
    int nc = 0;
    double p_state = 1.0 / num_states;

    for (int si = 0; si < num_states; si++) {
        constraints[nc].coeffs = calloc(n_vars, sizeof(double));
        for (int h = 0; h < num_dice; h++)
            constraints[nc].coeffs[si * num_dice + h] = 1.0;
        constraints[nc].sense = '=';
        constraints[nc].rhs = p_state;
        nc++;
    }

    for (int sig_i = 0; sig_i < num_signals; sig_i++) {
        int b_val = sig_arr[sig_i][1];
        int z_idx = num_x + sig_i;

        constraints[nc].coeffs = calloc(n_vars, sizeof(double));
        constraints[nc + 1].coeffs = calloc(n_vars, sizeof(double));
        constraints[nc].coeffs[z_idx] = -1.0;
        constraints[nc + 1].coeffs[z_idx] = -1.0;

        for (int si = 0; si < num_states; si++) {
            for (int h = 0; h < num_dice; h++) {
                int revealed[2];
                int ri = 0;
                for (int i = 0; i < num_dice; i++) {
                    if (i != h) revealed[ri++] = states[si][i];
                }
                if (revealed[0] > revealed[1]) {
                    int tmp = revealed[0]; revealed[0] = revealed[1]; revealed[1] = tmp;
                }
                if (revealed[0] == sig_arr[sig_i][0] && revealed[1] == sig_arr[sig_i][1]) {
                    int x_idx = si * num_dice + h;
                    constraints[nc].coeffs[x_idx] += b_val;
                    constraints[nc + 1].coeffs[x_idx] += states[si][h];
                }
            }
        }
        constraints[nc].sense = '<';
        constraints[nc].rhs = 0.0;
        constraints[nc + 1].sense = '<';
        constraints[nc + 1].rhs = 0.0;
        nc += 2;
    }

    double *objective = calloc(n_vars, sizeof(double));
    for (int si = 0; si < num_z; si++)
        objective[num_x + si] = 1.0;

    double result = solve_lp(n_vars, constraints, n_constraints, objective);

    for (int i = 0; i < n_constraints; i++) free(constraints[i].coeffs);
    free(constraints);
    free(objective);
    free(states);

    return result;
}

double p982_native(void) {
    return build_and_solve(3);
}
