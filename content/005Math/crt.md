# 中国剩余定理（Chinese Remainder Thoerem）

```c++
#include <cstdio>

const int M = 999911658;
const int m[] = {2, 3, 4679, 35617};
int a[4];

long long pow(long long a, long long n, long long p) {
    long long res = 1;
    for (; n; n >>= 1, a = a * a % p) if (n & 1) res = res * a % p;
    return res;
}

long long inv(long long n, long long p) {
    return pow(n, p - 2, p);
}

long long ChineseRemainderTheorem() {
    long long res = 0;
    for (int i = 0; i < 4; i++) {
        long long x = inv(M / m[i], m[i]) % M;
        long long temp = (x * (M / m[i]) % M + M) % M;
        (res += temp * a[i] % M) %= M;
    }
    return res;
}

int main() {
    for (int i = 0; i < 4; i++) scanf("%d", &a[i]);
    printf("%lld\n", ChineseRemainderTheorem());
    
    return 0;
}
```
