# Baby Step Giant Step

```c++
#include <cstdio>
#include <cmath>
#include <map>

long long pow(long long a, long long n, long long p) {
    long long res = 1;
    for (; n; n >>= 1, a = a * a % p) if (n & 1) res = res * a % p;
    return res;
}

long long inv(long long a, long long p) {
    return pow(a, p - 2, p);
}

long long bsgs(long long a, long long b, long long p) {
    a %= p, b %= p;
    if (a == 0) return b == 0 ? 1 : -1;

    std::map<long long, long long> map;

    long long m = std::ceil(std::sqrt(p)), t = 1;
    for (int i = 0; i < m; i++) {
        if (!map.count(t)) map[t] = i;
        t = t * a % p;
    }

    long long k = inv(t, p), w = b;
    for (int i = 0; i < m; i++) {
        if (map.count(w)) return i * m + map[w];
        w = w * k % p;
    }

    return -1;
}

int main () {
    long long a, b, p;
    scanf("%lld %lld %lld", &a, &b, &p);

    long long ans = bsgs(a, b, p);
    printf("%lld\n", ans);

    return 0;
}
```
