#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// 被積分関数
double func(double x, double v)
{
    double t;
    t = tanh(x / sqrt(v) + 1.0 / v);
    return t * t * exp(-x * x / 2.0) / sqrt(2.0 * M_PI);
}

// シンプソン法
double simpson_integration(double a, double b, int n, double v)
{
    if (n % 2 != 0) n++; 
    double h = (b - a) / n;
    double sum = func(a, v) + func(b, v);

    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        if (i % 2 == 0) {
            sum += 2.0 * func(x, v);
        } else {
            sum += 4.0 * func(x, v);
        }
    }
    return (h / 3.0) * sum;
}

// ★ここが消えていないかチェック！
int main()
{
    double v = 0.1;
    double a = 6.0;
    int n = 1000;

    double integral_val = simpson_integration(-a, a, n, v);
    double mse = 1.0 - integral_val;

    printf("v = %.1f のときの計算結果:\n", v);
    printf("積分値   : %f\n", integral_val);
    printf("MSE(0.1) : %f\n", mse);

    return 0;
}
