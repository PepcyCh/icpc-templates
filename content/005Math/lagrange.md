# 拉格朗日插值（Lagrange Interpolation Polynomial）

```c++
#include <cstdio>

const int MAXN = 1000005;
const int MOD = 1000000007;

long long qpow(long long a, long long n) {
    long long res = 1;
    for (; n; n >>= 1, a = a * a % MOD) if (n & 1) res = res * a % MOD;
    return res;
}

long long fact[MAXN], invFact[MAXN];
void init() {
    fact[0] = 1;
    for (int i = 1; i < MAXN; i++) fact[i] = fact[i - 1] * i % MOD;
    invFact[MAXN - 1] = qpow(fact[MAXN - 1], MOD - 2);
    for (int i = MAXN - 2; ~i; i--) invFact[i] = invFact[i + 1] * (i + 1) % MOD;
}

long long Lagrange(long long *a, int n, long long x) {
    if (x <= n) return a[x];

    static long long f1[MAXN], f2[MAXN];
    f1[0] = f2[n + 1] = 1;
    for (int i = 1; i <= n + 1; i++) f1[i] = f1[i - 1] * (x - i + 1) % MOD;
    for (int i = n; ~i; i--) f2[i] = f2[i + 1] * (x - i) % MOD;

    long long res = 0;
    for (int i = 0; i <= n; i++) {
        long long temp = a[i] * f1[i] % MOD * f2[i + 1] % MOD * invFact[i] % MOD * invFact[n - i] % MOD;
        (n - i) % 2 ? res -= temp : res += temp;
        res >= MOD ? res -= MOD : 0;
        res < 0 ? res += MOD : 0;
    }

    return res;
}

int main() {
    init();

    int n;
    scanf("%d", &n);

    static long long a[MAXN + 1];
    for (int i = 0; i <= n; i++) scanf("%lld", &a[i]);

    long long x;
    scanf("%lld", &x);
    long long ans = Lagrange(a, n, x);
    printf("%lld\n", ans);

    return 0;
}
```
