#include <stdio.h>
#include <math.h>

double f(double x)
{
    return x * x * x - x;
}

int main(void)
{
    double a = -0.5;
    double b = 0.5;
    double c;
    double true_root = 0.0;
    double error;
    int i;

    printf("方程式 f(x) = x^3 - x = 0\n");
    printf("求める解 : x2 = 0\n\n");

    printf("=== 二分法の初期値が満たすべき条件 ===\n");
    printf("1. a < b\n");
    printf("2. f(a) * f(b) < 0 （解を挟む区間であること）\n\n");

    printf("初期値\n");
    printf("a = %f\n", a);
    printf("b = %f\n", b);
    printf("f(a) = %f\n", f(a));
    printf("f(b) = %f\n", f(b));

    /* 初期条件の確認 */
    if (a >= b)
    {
        printf("\n条件 a < b を満たしていません。\n");
        return 1;
    }

    if (f(a) * f(b) >= 0)
    {
        printf("\n条件 f(a) * f(b) < 0 を満たしていません。\n");
        return 1;
    }

    printf("\n初期値は二分法の条件を満たしている。\n");

    printf("\n反復回数\t近似解\t\t絶対誤差\n");

    for (i = 0; i < 10; i++)
    {
        c = (a + b) / 2.0;
        error = fabs(true_root - c);

        printf("%d\t\t%.10f\t%.10f\n", i, c, error);

        /* 真の解を含む区間を選択 */
        if (f(a) * f(c) < 0)
        {
            b = c;
        }
        else
        {
            a = c;
        }
    }

    return 0;
}
