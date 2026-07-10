#include <stdio.h>
#include <math.h>

#define EPS 1.0e-15 

double Function(double x)
{
    return x * x * x - x;
}

void Bisection(double xl0, double xr0)
{
    int i, number;
    double x, xl, xr;
    double yl, yx; // 関数の値を保存する変数

    number = 10000;
    xl = xl0; 
    xr = xr0; 
    
    yl = Function(xl);

    // ループ条件: 最大回数まで、かつ区間の幅がEPSより大きい間
    for (i = 1; i <= number && (xr - xl) > EPS; i++) {
        x = (xl + xr) / 2.0;
        yx = Function(x); // 関数の計算は1回だけにシンプライズ

        // 真の解 0.0 との絶対誤差を表示
        printf("Iteration=%d, solution=%g, error=%g\n", i, x, fabs(x - 0.0)); 
        
        if (yx == 0.0) {
            break; // ぴったり解が見つかったら終了
        }

        // 符号が同じかどうかの判定（掛け算によるアンダーフローを防止）
        // signbit() は符号がマイナスの時に真を返す関数です
        if (signbit(yl) == signbit(yx)) {
            xl = x;
            yl = yx; // 左端が更新されたので、左端の関数値も更新
        } else {
            xr = x;
        }
    }
}

int main()
{
    // 初期条件: -1 < xl0 < 0 < xr0 < 1
    Bisection(-0.4, 0.6);
    return 0;
}
