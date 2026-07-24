/*
 * Active Learning 6
 * sigma^2 = 0.1, alpha = 1.8 のとき
 *   v = sigma^2 + alpha * MSE(v),   MSE(v) = v / (1 + v)
 * を満たす v > 0 を二分法で求める。
 *
 * 1. MSE(v) の関数形の情報を「二分法のロジック内」では使わず、
 *    f(v) = sigma^2 + alpha*MSE(v) - v の符号だけを使って解を求める。
 * 2. MSE(v) = v/(1+v) を代入して閉形式(2次方程式)を導出し、
 *    数値解と比較する。
 */

#include <stdio.h>
#include <math.h>

/* ---- 問題のパラメータ ---- */
#define SIGMA2 0.1
#define ALPHA  1.8

/* ---- MSE(v) の定義（ここにだけ関数形の情報がある） ---- */
double MSE(double v)
{
    return v / (1.0 + v);
}

/* ---- 方程式 v = sigma^2 + alpha*MSE(v) を f(v)=0 の形にしたもの ---- */
double f(double v)
{
    return SIGMA2 + ALPHA * MSE(v) - v;
}

/*
 * ---- 汎用二分法 ----
 * f の関数形（MSE(v)がどんな式か）には一切依存しない。
 * f(a)*f(b) < 0 という符号情報だけを使う。
 */
double bisection(double (*func)(double), double a, double b, double tol, int max_iter)
{
    double fa = func(a);
    double fb = func(b);
    double c = a, fc;

    if (fa * fb > 0.0) {
        fprintf(stderr, "Error: f(a) と f(b) の符号が同じです。区間を見直してください。\n");
        return NAN;
    }

    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        c = 0.5 * (a + b);
        fc = func(c);

        printf("iter=%2d  a=%.10f  b=%.10f  c=%.10f  f(c)=%.3e\n",
               iter, a, b, c, fc);

        if (fabs(fc) < tol || 0.5 * (b - a) < tol) {
            break;
        }

        if (fa * fc < 0.0) {
            b  = c;
            fb = fc;
        } else {
            a  = c;
            fa = fc;
        }
    }
    return c;
}

int main(void)
{
    double a = 0.0;      /* v > 0 の探索区間下限 */
    double b = 10.0;     /* 十分大きな上限（f(a)*f(b)<0 となるように） */
    double tol = 1e-10;
    int max_iter = 100;

    printf("=== 1. 二分法による数値解 ===\n");
    printf("f(v) = sigma^2 + alpha*MSE(v) - v,  sigma^2=%.3f, alpha=%.3f\n\n",
           SIGMA2, ALPHA);

    double v_num = bisection(f, a, b, tol, max_iter);

    printf("\n二分法による解: v = %.10f\n", v_num);

    /* ---- 2. MSE(v) = v/(1+v) を代入して閉形式を導出 ----
     *
     * v = sigma^2 + alpha * v/(1+v)
     * v(1+v) = sigma^2(1+v) + alpha*v
     * v + v^2 = sigma^2 + sigma^2*v + alpha*v
     * v^2 + (1 - sigma^2 - alpha) v - sigma^2 = 0
     *
     * これは A*v^2 + B*v + C = 0 の形（A=1）
     */
    double A = 1.0;
    double B = 1.0 - SIGMA2 - ALPHA;
    double C = -SIGMA2;

    double disc = B * B - 4.0 * A * C;
    double v1 = (-B + sqrt(disc)) / (2.0 * A);
    double v2 = (-B - sqrt(disc)) / (2.0 * A);

    printf("\n=== 2. 閉形式（解析解）との比較 ===\n");
    printf("v^2 + (%.4f) v + (%.4f) = 0\n", B, C);
    printf("判別式 = %.10f\n", disc);
    printf("解1: v = %.10f\n", v1);
    printf("解2: v = %.10f\n", v2);

    double v_analytic = (v1 > 0) ? v1 : v2; /* v > 0 の方を採用 */

    printf("\n解析解 (v > 0): v = %.10f\n", v_analytic);
    printf("数値解 (二分法): v = %.10f\n", v_num);
    printf("誤差 |v_num - v_analytic| = %.3e\n", fabs(v_num - v_analytic));

    return 0;
}
