# 自适应辛普森积分（Adaptive Simpson's Method）

```c++
#include <cstdio>
#include <cmath>

double f(double x) {
    return std::sqrt(x);
}

double simpson(double l, double r) {
    double mid = (l + r) / 2.0;
    return (f(l) + 4.0 * f(mid) + f(r)) * (r - l) / 6.0;
}

double integral(double l, double r, double eps) {
    double mid = (l + r) / 2.0;
    double ST = simpson(l, r), SL = simpson(l, mid), SR = simpson(mid, r);
    if (std::abs(SL + SR - ST) <= 15.0 * eps) return SL + SR + (SL + SR - ST) / 15.0;
    return integral(l, mid, eps / 2.0) + integral(mid, r, eps / 2.0);
}

int main() {
    double l, r, eps;
    scanf("%lf %lf %lf", &l, &r, &eps);
    printf("%.6f\n", integral(l, r, eps));
    
    return 0;
}
```
