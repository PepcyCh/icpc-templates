# 组合数与 Lucas 定理

```c++
#include <cstdio>

const int MOD = 1000003;

int fact[MOD], inv[MOD];

long long pow(long long a, long long n) {
    long long res = 1;
    for (; n; n >>= 1, a = a * a % MOD) if (n & 1) res = res * a % MOD;
    return res;
}

void prepare() {
    fact[0] = 1;
    for (int i = 1; i < MOD; i++) fact[i] = (long long) fact[i - 1] * i % MOD;

    inv[MOD - 1] = pow(fact[MOD - 1], MOD - 2);
    for (int i = MOD - 2; i; i--) inv[i] = (long long) inv[i + 1] * (i + 1) % MOD;
}

int combi(int n, int m) {
    if (m > n) return 0;
    if (n < MOD) return (long long) fact[n] * inv[m] % MOD * inv[n - m] % MOD;
    return (long long) combi(n / MOD, m / MOD) * combi(n % MOD, m % MOD) % MOD;
}

int main() {
    prepare();

    return 0;
}
```
