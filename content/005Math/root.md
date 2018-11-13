# 原根

```c++
#include <cstdio>

const int MAXN = 1000005;

int prime[MAXN], primeCnt;

void sieve() {
    static bool notPrime[MAXN];
    notPrime[0] = notPrime[1] = true;
    primeCnt = 0;

    for (int i = 2; i < MAXN; i++) {
        if (!notPrime[i]) prime[primeCnt++] = i;

        for (int j = 0; j < primeCnt && i * prime[j] < MAXN; j++)
            notPrime[i * prime[j]] = true;
    }
}

int pow(int a, int n, int p) {
    int res = 1;
    for (; n; n >>= 1, a = a * a % p) if (n & 1) res = res * a % p;
    return res;
}

int getRoot(int p) {
    for (int g = 2, pp = p - 1; ; g++) {
        bool flag = true;
        for (int i = 0; i < primeCnt && prime[i] < p; i++) {
            if (pp % prime[i] == 0 && pow(g, pp / prime[i], p) == 1) {
                flag = false;
                break;
            }
        }

        if (flag) return g;
    }
}

int main() {

    return 0;
}
```
