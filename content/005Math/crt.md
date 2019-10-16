# 线性同余方程

```c++
#include <cstdio>

const int MAXN = 100005;

void exgcd(long long a, long long b, long long &g, long long &x, long long &y) {
    if (b == 0) x = 1, y = 0, g = a;
    else exgcd(b, a % b, g, y, x), y -= a / b * x;
}

int mod[MAXN], rem[MAXN];
long long solveCongruence(int n) {
    long long res = 0, K = 1;
    for (int i = 0; i < n; i++) {
        long long x, y, g;
        exgcd(K, mod[i], g, x, y);
        if ((rem[i] - res) % g) return -1;
        x = (x * (rem[i] - res) / g + mod[i] / g) % (mod[i] / g);
        y = (K / g *mod[i]);
        res = ((x * K + res) % y + y) % y;
        K = y;
    }
    return res;
}

int main() {
    int n;
    scanf("%d", &n);
    for (int i = 0; i < n; i++) scanf("%d %d", &mod[i], &rem[i]);
    long long ans = solveCongruence(n);
    printf("%lld\n", ans);
    
    return 0;
}
```
