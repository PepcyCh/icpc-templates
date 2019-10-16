# 欧拉筛（Euler Seive）

```c++
#include <cstdio>

const int MAXN = 10000005;

int prime[MAXN], primeCnt, mu[MAXN], phi[MAXN], d[MAXN];
long long s[MAXN];
bool notPrime[MAXN];

void sieve() {
    static int minFact[MAXN], minPow[MAXN];

    notPrime[0] = notPrime[1] = true;
    mu[1] = phi[1] = d[1] = s[1] = 1;
    minFact[1] = minPow[1] = 1;

    for (int i = 2; i < MAXN; i++) {
        if (!notPrime[i]) {
            prime[primeCnt++] = i;
            mu[i] = -1;
            phi[i] = i - 1;
            d[i] = 2;
            s[i] = i + 1;
            minFact[i] = i;
            minPow[i] = 1;
        }

        for (int j = 0; j < primeCnt && i * prime[j] < MAXN; j++) {
            notPrime[i * prime[j]] = true;
            if (i % prime[j] == 0) {
                mu[i * prime[j]] = 0;
                phi[i * prime[j]] = phi[i] * prime[j];
                d[i * prime[j]] = d[i] / (minPow[i] + 1) * (minPow[i] + 2);
                if (i == minFact[i]) {
                    s[i * prime[j]] = s[i] + i * prime[j];
                } else {
                    s[i * prime[j]] = s[i / minFact[i]] * s[prime[j] * minFact[i]];
                }
                minPow[i * prime[j]] = minPow[i] + 1;
                minFact[i * prime[j]] = minFact[i] * prime[j];
                break;
            }
            mu[i * prime[j]] = -mu[i];
            phi[i * prime[j]] = phi[i] * (prime[j] - 1);
            d[i * prime[j]] = d[i] * 2;
            s[i * prime[j]] = s[i] * s[prime[j]];
            minPow[i * prime[j]] = 1;
            minFact[i * prime[j]] = prime[j];
        }
    }
}

int main() {
    sieve();

    return 0;
}
```
