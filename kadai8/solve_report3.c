#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SIGMA2 0.08
#define X_MIN -8.0
#define X_MAX 8.0
#define SIMPSON_N 4000
#define V_SCAN_STEP 0.001
#define BISECTION_TOL 1.0e-10
#define BISECTION_MAX_ITER 100

static double X[SIMPSON_N + 1];
static double WEIGHTED_DENSITY[SIMPSON_N + 1];
static double SIMPSON_FACTOR;

static void make_simpson_data(void) {
    if (SIMPSON_N <= 0 || SIMPSON_N % 2 != 0) {
        fprintf(stderr, "SIMPSON_N must be a positive even integer\n");
        exit(EXIT_FAILURE);
    }

    double dx = (X_MAX - X_MIN) / SIMPSON_N;
    double inv_sqrt_2pi = 1.0 / sqrt(2.0 * M_PI);

    for (int i = 0; i <= SIMPSON_N; ++i) {
        double x = X_MIN + i * dx;
        X[i] = x;
        double weight = 1.0;
        if (i > 0 && i < SIMPSON_N) {
            if (i % 2 == 1) {
                weight = 4.0;
            } else {
                weight = 2.0;
            }
        }
        double normal_density = exp(-0.5 * x * x) * inv_sqrt_2pi;
        WEIGHTED_DENSITY[i] = weight * normal_density;
    }
    SIMPSON_FACTOR = dx / 3.0;
}

static double mse(double v) {
    if (v <= 0.0) {
        fprintf(stderr, "v must be positive\n");
        exit(EXIT_FAILURE);
    }

    double sum = 0.0;
    double inv_sqrt_v = 1.0 / sqrt(v);
    double inv_v = 1.0 / v;

    for (int i = 0; i <= SIMPSON_N; ++i) {
        double z = X[i] * inv_sqrt_v + inv_v;
        double t = tanh(z);
        sum += WEIGHTED_DENSITY[i] * (t * t);
    }
    return 1.0 - SIMPSON_FACTOR * sum;
}

static double equation(double v, double alpha) {
    return v - SIGMA2 - alpha * mse(v);
}

static double bisect(double alpha, double left, double right) {
    double f_left = equation(left, alpha);
    double f_right = equation(right, alpha);

    if (fabs(f_left) < BISECTION_TOL) {
        return left;
    }
    if (fabs(f_right) < BISECTION_TOL) {
        return right;
    }
    if (f_left * f_right > 0.0) {
        fprintf(stderr, "The interval does not bracket a root: [%.12f, %.12f]\n", left, right);
        exit(EXIT_FAILURE);
    }

    for (int iter = 0; iter < BISECTION_MAX_ITER; ++iter) {
        double middle = 0.5 * (left + right);
        double f_middle = equation(middle, alpha);

        if (fabs(f_middle) < BISECTION_TOL || right - left < BISECTION_TOL) {
            return middle;
        }
        if (f_left * f_middle <= 0.0) {
            right = middle;
            f_right = f_middle;
        } else {
            left = middle;
            f_left = f_middle;
        }
    }
    return 0.5 * (left + right);
}

static int find_all_roots(double alpha,
                          const double *v_grid,
                          const double *mse_grid,
                          int grid_size,
                          double *roots,
                          int max_roots) {
    int last = 0;
    double v_max = SIGMA2 + alpha;
    for (int i = 0; i < grid_size; ++i) {
        if (v_grid[i] > v_max) {
            last = i;
            break;
        }
    }
    if (last == 0) {
        last = grid_size;
    }

    int root_count = 0;
    double f_prev = v_grid[0] - SIGMA2 - alpha * mse_grid[0];

    for (int i = 0; i < last - 1; ++i) {
        double f_curr = v_grid[i + 1] - SIGMA2 - alpha * mse_grid[i + 1];

        if (f_prev == 0.0 || f_curr == 0.0 || f_prev * f_curr < 0.0) {
            double root;
            if (fabs(f_prev) < 1.0e-13) {
                root = v_grid[i];
            } else if (fabs(f_curr) < 1.0e-13) {
                root = v_grid[i + 1];
            } else {
                root = bisect(alpha, v_grid[i], v_grid[i + 1]);
            }
            if (root_count == 0 || fabs(root - roots[root_count - 1]) > 10.0 * BISECTION_TOL) {
                if (root_count < max_roots) {
                    roots[root_count++] = root;
                }
            }
        }
        f_prev = f_curr;
    }
    return root_count;
}

static void print_three_root_ranges(const double *alphas,
                                    const int *counts,
                                    int n_alphas) {
    int start_index = -1;
    int previous_index = -1;
    int found = 0;

    for (int i = 0; i < n_alphas; ++i) {
        if (counts[i] == 3) {
            if (start_index < 0) {
                start_index = i;
            } else if (i - previous_index != 1) {
                printf("%.2f <= alpha <= %.2f\n", alphas[start_index], alphas[previous_index]);
                start_index = i;
            }
            previous_index = i;
            found = 1;
        }
    }
    if (found && start_index >= 0) {
        printf("%.2f <= alpha <= %.2f\n", alphas[start_index], alphas[previous_index]);
    }
    if (!found) {
        printf("該当範囲なし\n");
    }
}

int main(void) {
    make_simpson_data();

    const int alpha_count = 901;
    double *alphas = malloc(sizeof(double) * alpha_count);
    if (!alphas) {
        fprintf(stderr, "Memory allocation failed\n");
        return EXIT_FAILURE;
    }
    for (int i = 0; i < alpha_count; ++i) {
        alphas[i] = 1.0 + i * 0.01;
    }

    double v_max = SIGMA2 + 10.0;
    int scan_count = (int)round((v_max - SIGMA2) / V_SCAN_STEP);
    int grid_size = scan_count + 1;

    double *v_grid = malloc(sizeof(double) * grid_size);
    double *mse_grid = malloc(sizeof(double) * grid_size);
    int *root_counts = malloc(sizeof(int) * alpha_count);
    double *roots_buffer = malloc(sizeof(double) * alpha_count * 10);

    if (!v_grid || !mse_grid || !root_counts || !roots_buffer) {
        fprintf(stderr, "Memory allocation failed\n");
        return EXIT_FAILURE;
    }

    for (int i = 0; i < grid_size; ++i) {
        v_grid[i] = SIGMA2 + i * V_SCAN_STEP;
    }

    printf("探索用の MSE(v) を計算中...\n");
    for (int i = 0; i < grid_size; ++i) {
        mse_grid[i] = mse(v_grid[i]);
    }

    FILE *file = fopen("solutions.csv", "w");
    if (!file) {
        fprintf(stderr, "Failed to open solutions.csv for writing\n");
        return EXIT_FAILURE;
    }
    fprintf(file, "alpha,root_number,v\n");

    for (int i = 0; i < alpha_count; ++i) {
        double alpha = alphas[i];
        double *roots = roots_buffer + i * 10;
        int n_roots = find_all_roots(alpha, v_grid, mse_grid, grid_size, roots, 10);
        root_counts[i] = n_roots;
        for (int j = 0; j < n_roots; ++j) {
            fprintf(file, "%.2f,%d,%.10f\n", alpha, j + 1, roots[j]);
        }
    }
    fclose(file);

    printf("\nalpha ごとの正の解\n");
    for (int i = 0; i < alpha_count; ++i) {
        double alpha = alphas[i];
        double *roots = roots_buffer + i * 10;
        printf("alpha=%5.2f: %d solution(s)", alpha, root_counts[i]);
        if (root_counts[i] > 0) {
            printf(": ");
            for (int j = 0; j < root_counts[i]; ++j) {
                printf("%.10f", roots[j]);
                if (j < root_counts[i] - 1) {
                    printf(", ");
                }
            }
        }
        printf("\n");
    }

    printf("\n正の解が3個存在する alpha の範囲（刻み幅 0.01）\n");
    print_three_root_ranges(alphas, root_counts, alpha_count);

    printf("\nCSV: solutions.csv\n");
    printf("(plot は C 版には含まれていません)\n");

    free(alphas);
    free(v_grid);
    free(mse_grid);
    free(root_counts);
    free(roots_buffer);
    return EXIT_SUCCESS;
}
