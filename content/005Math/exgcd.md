# 欧几里得算法（Euclidean Algorithm）

```c++
#include <cstdio>

int gcd(int a, int b) {
    return b ? gcd(b, a % b) : a;
}

void exgcd(int a, int b, int &g, int &x, int &y) {
    if (b == 0) x = 1, y = 0, g = a;
    else exgcd(b, a % b, g, y, x), y -= a / b * x;
}

int main() {

    return 0;
}
```
