#include <stdio.h>
#include <math.h>

// M_PIを使うための定義（環境によって未定義の場合の対策）
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// 被積分関数の定義（スライド3ページ目より）
double func(double x, double v) {
    double t = tanh(x / sqrt(v) + 1.0 / v);
    return t * t * exp(-x * x / 2.0) / sqrt(2.0 * M_PI);
}

// シンプソン法を使用して数値積分を計算する関数
double SimpInt(int n, double a, double b, double v) {
    int i;
    double S, h;
    S = 0.0;
    h = (b - a) / (double)n;
    for (i = 0; i < n; i++) {
        S += func(i * h + a, v) + 4.0 * func((i + 0.5) * h + a, v) + func((i + 1.0) * h + a, v);
    }
    S *= h / 6.0;
    return S;
}

// 誤差評価のための関数（スライド4ページ目の上界式）
// 2 * ∫_{a}^{∞} (1/√(2π)) * exp(-x^2/2) dx = 1 - erf(a/√2)
double error_estimation(double a) {
    return 1.0 - erf(a / sqrt(2.0));
}

int main() {
    double v = 0.1;
    double a = 0.0; // 初期の積分区間の上限（0から探索開始）
    int n = 5000;  // 分割数
    double tolerance = 1e-7; // 許容誤差（スライドより 10^-7）

    // 誤差が10^-7以下になる最小のa>0を0.1刻みで探索
    while (error_estimation(a) > tolerance) {
        a += 0.1;
    }

    // 数値積分を計算（積分区間は [-a, a]）
    double integral = SimpInt(n, -a, a, v);

    // MSE(v)を計算
    double mse = 1.0 - integral;

    // 結果を表示
    printf("MSE(%.1f) = %.10f\n", v, mse);
    printf("Integral range: [-%.1f, %.1f]\n", a, a);
    printf("Estimated error: %e\n", error_estimation(a));
    printf("Final value of a: %.1f\n", a);

    return 0;
}
