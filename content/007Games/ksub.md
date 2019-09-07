# K 倍动态减法

一堆石子有 n 个，第一次可以取 1 ~ n - 1 个石子，之后至多可取 k 倍上次所取的石子数。

k 时的先手必败集合是一个 k 次的递推。

k = 1: a_n = 2 a_{n - 1}

k = 2: a_n = a_{n - 1} + a_{n - 2}

k = 3: a_n = a_{n - 1} + a_{n - 2} - a_{n - 3} + 1

```c++
#include <cstdio>

const int MAXN = 2000005;

int a[MAXN], b[MAXN];

int main() {
    int n, k; // n <= 100,000,000, k <= 100,000
    scanf("%d %d", &n, &k);

    a[0] = b[0] = 1;
    int i = 0, j = 0;
    while (a[i] < n) {
        ++i;
        a[i] = b[i - 1] + 1;
        while (a[j + 1] * k < a[i]) ++j;
        if (a[j] * k < a[i]) b[i] = a[i] + b[j];
        else b[i] = a[i];
    }

    if (a[i] == n) puts("Lose");
    else {
        int ans = 0;
        while (n) {
            if(n >= a[i]) {
                n -= a[i];
                ans = a[i];
            }
            --i;
        }
        printf("%d\n",ans);
    }

    return 0;
}
```
