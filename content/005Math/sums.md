# 杜教筛

```c++
#include <cstdio>
#include <map>

const int MAXNN = 1600000;

long long phi[MAXNN];
int prime[MAXNN], primeCnt;
bool notPrime[MAXNN];

void seive() {
    phi[1] = 1;
    notPrime[0] = notPrime[1] = true;
    for (int i = 2; i < MAXNN; i++) {
        if (!notPrime[i]) {
            prime[++primeCnt] = i;
            phi[i] = i - 1;
        }

        for (int j = 1; j <= primeCnt && i * prime[j] < MAXNN; j++) {
            notPrime[i * prime[j]] = true;
            if (i % prime[j] == 0) {
                phi[i * prime[j]] = phi[i] * prime[j];
                break;
            } else phi[i * prime[j]] = phi[i] * (prime[j] - 1);
        }
    }

    for (int i = 2; i < MAXNN; i++) phi[i] += phi[i - 1];
}

/*
 *  h = f * g
 *  g(1) F(n) = H(n) - \sum_{i = 2}^{n} g(i) F(\lfloor \frac{n}{i} \rfloor)
 */
long long pSumPhi(int n) {
    if (n < MAXNN) return ::phi[n];

    static std::map<int, long long> phi;
    std::map<int, long long>::iterator it;

    if ((it = phi.find(n)) != phi.end()) return it->second;

    long long res = (long long) n * (n + 1) / 2;
    for (int i = 2, last; i <= n; i = last + 1) {
        last = n / (n / i);
        res -= (last - i + 1) * pSumPhi(n / i);
    }
    return phi[n] = res;
}

int main() {
    seive();

    int n;
    scanf("%d", &n);
    printf("%lld\n", pSumPhi(n));

    return 0;
}
```
