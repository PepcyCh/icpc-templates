# Min25 筛

```c++
#include <cstdio>
#include <cmath>

const int MAXN = 100005; // sqrt n
const int MOD = 1000000007;
int INV2, INV6;

long long qpow(long long a, long long n) {
    long long res = 1;
    for (; n; n >>= 1, a = a * a % MOD) if (n & 1) res = res * a % MOD;
    return res;
}

int prime[MAXN], primeCnt, lpf[MAXN];
long long sprime[MAXN], ssprime[MAXN];
bool notPrime[MAXN];

void sieve() {
    for (int i = 2; i < MAXN; i++) {
        if (!notPrime[i]) {
            prime[++primeCnt] = i; // index from 1
            sprime[primeCnt] = (sprime[primeCnt - 1] + i) % MOD;
            ssprime[primeCnt] = (ssprime[primeCnt - 1] + 1ll * i * i) % MOD;
            lpf[i] = primeCnt;
        }
        for (int j = 1; j <= primeCnt && i * prime[j] < MAXN; j++) {
            notPrime[i * prime[j]] = true;
            if (i % prime[j] == 0) break;
            lpf[i * prime[j]] = j;
        }
    }
}

class Min25 {
private:
    long long G[MAXN][3], Fprime[MAXN];
    int list[MAXN], cnt;

    int le[MAXN], ge[MAXN], lim, n;
    inline int &id(int v) { return v <= lim ? le[v] : ge[n / v]; }

    void init(int n) {
        cnt = 0;
        for (int i = 1, j, v; i <= n; i = n / j + 1) {
            j = n / i;
            v = j % MOD;
            list[++cnt] = j;
            id(j) = cnt;
            // sum from 2 to j (assuming all the values are calculated as if they are prime)
            // i^0
            G[cnt][0] = v ? v - 1 : MOD - 1;
            // i^1
            G[cnt][1] = (2ll + v) * (v - 1ll) % MOD * INV2 % MOD;
            // i^2
            G[cnt][2] = (1ll * v * (v + 1) % MOD * (2 * v + 1) % MOD * INV6 % MOD - 1 + MOD) % MOD;
        }
    }

    void calcFprime() {
        for (int k = 1; k <= primeCnt && prime[k] <= lim; k++) {
            int p = prime[k];
            long long sqrp = 1ll * p * p;
            for (int i = 1; list[i] >= sqrp; i++) {
                int v = list[i] / p;
                int id = this->id(v);
                G[i][0] = (G[i][0] - G[id][0] + k - 1 + MOD) % MOD;
                G[i][1] = (G[i][1] - 1ll * p * (G[id][1] - sprime[k - 1] + MOD) % MOD + MOD) % MOD;
                G[i][2] = (G[i][2] - 1ll * p * p % MOD * (G[id][2] - ssprime[k - 1] + MOD) % MOD + MOD) % MOD;
            }
        }
        for (int i = 1; i <= cnt; i++) {
            // prefix sum of values at primes
            Fprime[i] = (G[i][1] - G[i][0] + MOD) % MOD;
        }
    }

    int fp(int p, int c) {
        // value at power of prime
        return qpow(p, c - 1) * (p - 1) % MOD;
    }

    long long F(int K, int n) {
        if (n < prime[K] || n <= 1) return 0;
        int id = this->id(n);
        long long ans = Fprime[id] - (sprime[K - 1] - (K - 1));
        ans = (ans % MOD + MOD) % MOD;
        for (int i = K; i <= primeCnt && 1ll * prime[i] * prime[i] <= n; i++) {
            long long pw = prime[i], pw2 = pw * pw;
            for (int c = 1; pw2 <= n; pw = pw2, pw2 *= prime[i], c++) {
                ans = (ans + fp(prime[i], c) * F(i + 1, n / pw) + fp(prime[i], c + 1)) % MOD;
            }
        }
        return ans;
    }

public:
    long long solve(int n) {
        this->n = n;
        lim = sqrt(n);
        init(n);
        calcFprime();
        long long res = F(1, n) + 1; // remember to add f(1)
        return res >= MOD ? res -= MOD : res;
    }
} min25;

int main() {
    sieve();
    INV2 = qpow(2, MOD - 2);
    INV6 = qpow(6, MOD - 2);

    int n;
    scanf("%d", &n);
    long long ans = min25.solve(n);
    printf("%lld\n", ans);

    return 0;
}
```
