# 多项式 GCD

```c++
#include <vector>

const int MOD = 1000000007;

long long qpow(long long a, long long n) {
    long long res = 1;
    for (; n; n >>= 1, a = a * a % MOD) if (n & 1) res = res * a % MOD;
    return res;
}

// a[i] is the coefficient of x^i
std::vector<int> polyGcd(std::vector<int> a, const std::vector<int> &b) {
    if (b.empty()) return a;
    int t = a.size() - b.size();
    for (int i = 0; i <= t; i++) {
        long long temp = a[a.size() - 1 - i] * qpow(b[b.size() - 1], MOD - 2) % MOD;
        for (int j = 0; j < b.size(); j++)
            a[a.size() - 1 - i - j] = (a[a.size() - 1 - i - j] - temp * b[b.size() - 1 - j] % MOD + MOD) % MOD;
    }
    int p = -1;
    for (int i = a.size() - 1; ~i; i--) if (a[i]) {
        p = i;
        break;
    }
    a.resize(p + 1);
    return polyGcd(b, a);
}

int main() {

    return 0;
}
```
