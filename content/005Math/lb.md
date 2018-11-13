# 线性基（Linear Base）

### 在线版本

```c++
#include <cstdio>

const int MAXL = 50;

struct LinearBasis {
    long long a[MAXL + 1];

    bool extend(long long t) {
        for (int i = MAXL; ~i; i--) {
            if (!(t & (1ll << i))) continue;
            if (!t) return false;

            if (a[i]) t ^= a[i];
            else {
                for (int j = 0; j < i; j++) if (t & (1ll << j)) t ^= a[j];
                for (int j = i + 1; j <= MAXL; j++) if (a[j] & (1ll << i)) a[j] ^= t;
                a[i] = t;
                return true;
            }
        }

        return false;
    }

    long long queryMax() {
        long long res = 0;
        for (int i = 0; i <= MAXL; i++) res ^= a[i];
        return res;
    }
} lb;

int main() {
    int n;
    scanf("%d", &n);

    for (int i = 0; i < n; i++) {
        long long x;
        scanf("%lld", &x);
        lb.extend(x);
    }

    printf("%lld\n", lb.queryMax());

    return 0;
}
```

### 离线版本

```c++
#include <cstdio>
#include <vector>
#include <algorithm>

const int MAXN = 100005;
const int MAXL = 50;

struct LinearBasis {
    std::vector<long long> v;
    int n;

    void init(long long *x, int n) {
        this->n = n;
        static long long a[MAXL + 1];

        for (int i = 0; i < n; i++) {
            long long t = x[i];

            for (int j = MAXL; ~j; j--) {
                if (!(t & (1ll << j))) continue;
                if (!t) break;

                if (a[j]) t ^= a[j];
                else {
                    for (int k = 0; k < j; k++) if (t & (1ll << k)) t ^= a[k];
                    for (int k = j + 1; k <= MAXL; k++) if (a[k] & (1ll << j)) a[k] ^= t;
                    a[j] = t;
                    break;
                }
            }
        }

        v.clear();
        for (int i = 0; i <= MAXL; i++) if (a[i]) v.push_back(a[i]);
    }

    long long query(int k) {
        if (v.size() != n) k--;
        if (k >= (1ll << v.size())) return -1;

        long long res = 0;
        for (int i = 0; i < v.size(); i++) if (k & (1ll << i)) res ^= v[i];
        return res;
    }
} lb;

int main() {
    int n;
    scanf("%d", &n);

    static long long a[MAXN];
    for (int i = 0; i < n; i++) scanf("%lld", &a[i]);
    lb.init(a, n);

    int q;
    scanf("%d", &q);
    while (q--) {
        long long k;
        scanf("%lld", &k);
        printf("%lld\n", lb.query(k));
    }
    return 0;
}
```
