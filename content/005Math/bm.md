# Berlekamp-Massey

最短递推式。

```c++
#include <cstdio>
#include <vector>
#include <algorithm>

const int MAXN = 2505;
const int MOD = 998244353;

long long qpow(long long a, long long n) {
    long long res = 1;
    for (; n; n >>= 1, a = a * a % MOD) if (n & 1) res = res * a % MOD;
    return res;
}

class BerlekampMassey {
public:
    std::vector<int> solve(int *x, int n) {
        int L = 0, m = 1, b = 1;
        C = B = { 1 };
        for (int i = 0; i < n; i++) {
            int d = x[i];
            for (int j = 1; j <= L; j++) d = (d + 1ll * C[j] * x[i - j]) % MOD;
            if (d == 0) {
                ++m;
            } else if (2 * L <= i) {
                T = C;
                C.resize(i + 2 - L);
                int ib = qpow(b, MOD - 2);
                d = MOD - d;
                for (int i = 0; i < B.size(); i++) C[i + m] = (C[i + m] + 1ll * d * ib % MOD * B[i]) % MOD;
                L = i + 1 - L;
                B.swap(T);
                b = MOD - d;
                m = 1;
            } else {
                int ib = qpow(b, MOD - 2);
                d = MOD - d;
                for (int i = 0; i < B.size(); i++) C[i + m] = (C[i + m] + 1ll * d * ib % MOD * B[i]) % MOD;
                ++m;
            }
        }
        for (int i = 1; i <= L; i++) C[i] = (MOD - C[i]) % MOD;
        return C;
    }

private:
    std::vector<int> C, B, T;
} bm;

int x[MAXN];

int main() {
    int n;
    scanf("%d", &n);
    for (int i = 0; i < n; i++) scanf("%d", &x[i]);

    auto f = bm.solve(x, n);
    for (int i = 1; i < f.size(); i++) { // from 1
        printf("%d%c", f[i], " \n"[i == f.size() - 1]);
    }
    
    return 0;
}
```
