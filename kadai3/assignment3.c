#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// 被積分関数: f(x) = 1 / (1 + x)
// これを 0 から 1 まで積分すると log(2) になります
double func(double x)
{
    return 1.0 / (1.0 + x);
}

// シンプソン則による積分計算
double SimpInt(int n, double a, double b)
{
    int i;
    double S, h;

    S = 0.0; 
    h = (b - a) / (double)n;
    for (i = 0; i < n; i++) {
        S += func(i * h + a) + 4.0 * func((i + 0.5) * h + a) + func((i + 1.0) * h + a);
    }
    S *= h / 6.0;

    return S;
}

int main()
{
    int n;
    double a, b, S, S0, err;

    S0 = log(2.0); 
    a = 0.0; 
    b = 1.0;
    
    for (n = 10; n <= 1e+7; n *= 10) {
        S = SimpInt(n, a, b); 
        err = fabs((S - S0) / S0);
        printf("n=%8d Simpson=%.15f error=%g\n", n, S, err);
    }

    return 0;
}
