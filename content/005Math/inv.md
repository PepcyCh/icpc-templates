# 线性预处理逆元

```c++
#include <cstdio>

const int MAXN = 3000005;

int fact[MAXN], invFact[MAXN], inv[MAXN];

void exgcd(int a, int b, int &x, int &y) {
    if (b == 0) x = 1, y = 0;
    else exgcd(b, a % b, y, x), y -= x * (a / b);
}

int getInv(int x, int mod) {
    int res, temp;
    exgcd(x, mod, res, temp);
    return res;
}

int main() {
    int n, p;
    scanf("%d %d", &n, &p);

    inv[1] = 1;
    for (int i = 2; i <= n; i++) inv[i] = (long long) (p - p / i) * inv[p % i] % p;

    fact[0] = 1;
    for (int i = 1; i <= n; i++) fact[i] = (long long) fact[i - 1] * i % p;
    invFact[n] = getInv(fact[n], p);
    for (int i = n - 1; i; i--) invFact[i] = (long long) invFact[i + 1] * (i + 1)% p;

    return 0;
}
```
