# 常系数线性齐次递推

$$
h_k = \sum_{i = 1}^{k} a_{i - 1}h_{k - i}
$$

```c++
#include <cstdio>
#include <algorithm>

const int MOD = 1000000007;
const int MAXK = 2005;

long long a[MAXK];
void modMul(long long *A, long long *B, int n, long long *r) {
    static long long res[MAXK << 1];
    std::fill(res, res + (n << 1), 0);

    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
        res[i + j] += A[i] * B[j] % MOD;
        res[i + j] >= MOD ? res[i + j] -= MOD : 0;
    }

    for (int i = (n << 1) - 2; i >= n; i--) if (res[i]) for (int j = n - 1; ~j; j--) {
        res[i - n + j] += res[i] * a[j] % MOD;
        res[i - n + j] >= MOD ? res[i - n + j] -= MOD : 0;
    }

    for (int i = 0; i < n; i++) r[i] = res[i];
}

int main() {
    int n, k;
    scanf("%d %d", &n, &k);

    static long long h[MAXK];
    for (int i = 0; i < k; i++) {
        scanf("%lld", &a[i]);
        a[i] < 0 ? a[i] += MOD : 0;
    }
    std::reverse(a, a + k);
    a[k] = 1;

    for (int i = 0; i < k; i++) {
        scanf("%lld", &h[i]);
        h[i] < 0 ? h[i] += MOD : 0;
    }

    if (n < k) return printf("%lld\n", h[n]), 0;

    static long long m[MAXK], t[MAXK];
    m[0] = t[1] = 1;

    for (int i = n; i; i >>= 1, modMul(t, t, k, t)) if (i & 1) modMul(m, t, k, m);

    long long hn = 0;
    for (int i = 0; i < k; i++) hn = (hn + m[i] * h[i] % MOD) % MOD;

    printf("%lld\n", hn);

    return 0;
}
```
