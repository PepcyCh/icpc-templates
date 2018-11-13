# 高斯-约当消元（Gauss-Jordan Elimination）

```c++
#include <cstdio>
#include <cmath>
#include <algorithm>

const int MAXN = 105;
const double EPS = 1e-8;

int dcmp(double a, double b = 0.0) {
    double d = a - b;
    return std::abs(d) <= EPS ? 0 : (d > 0 ? 1 : -1);
}

namespace GuassJordan {
    double a[MAXN][MAXN];

    bool solve(int n) {
        for (int i = 0; i < n; i++) {
            int max = i;
            for (int j = i + 1; j < n; j++)
                if (dcmp(std::abs(a[j][i]), std::abs(a[max][i])) > 0) max = j;

            if (!dcmp(a[max][i])) return false;
            if (max != i) for (int j = i; j <= n; j++) std::swap(a[i][j], a[max][j]);

            for (int j = 0; j < n; j++) if (i != j) for (int k = n; k >= i; k--)
                a[j][k] -= a[i][k] / a[i][i] * a[j][i];
        }

        return true;
    }
}

int main() {
    int n;
    scanf("%d", &n);

    using GuassJordan::a;

    for (int i = 0; i < n; i++) for (int j = 0; j <= n; j++)
        scanf("%lf", &a[i][j]);

    if (!GuassJordan::solve(n)) puts("-1");
    else for (int i = 0; i < n; i++) printf("%.4lf\n", a[i][n] / a[i][i]);

    return 0;
}
```
