# 二次剩余（Quadratic Residue）

```c++
#include <cstdio>
#include <tuple>

long long qpow(long long a, long long n, long long p) {
    long long res = 1;
    for (; n; n >>= 1, a = a * a % p) if (n & 1) res = res * a % p;
    return res;
}

int modSqrt(int a, int p) {
    if (p == 2) return a % p;
    int x;
    if (qpow(a, (p - 1) >> 1, p) == 1) {
        if (p % 4 == 3) {
            x = qpow(a, (p + 1) >> 2, p);
        } else {
            long long w;
            for (w = 1; qpow((w * w - a + p) % p, (p - 1) >> 1, p) == 1; w++) {}
            long long b0 = w, b1 = 1;
            w = (w * w - a + p) % p;
            long long r0 = 1, r1 = 0;
            int exp = (p + 1) >> 1;
            for (; exp; exp >>= 1) {
                if (exp & 1)
                    std::tie(r0, r1) = std::make_tuple((r0 * b0 + r1 * b1 % p * w) % p, (r0 * b1 + r1 * b0) % p);
                std::tie(b0, b1) = std::make_tuple((b0 * b0 + b1 * b1 % p * w) % p, 2 * b0 * b1 % p);
            }
            x = r0;
        }
        if (x * 2 > p) x = p - x;
        return x;
    }
    return -1;
}

int main() {
    int n, p;
    scanf("%d %d", &n, &p);
    int ans = modSqrt(n, p);
    printf("%d\n", ans);

    return 0;
}
```
