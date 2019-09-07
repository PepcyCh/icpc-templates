# Nim-3

三人版的 Nim 取石子游戏，取完的人为一位，他的下一人为三位。

先手三位的条件是每个二进制位的 1 的个数为 3 的倍数。

```c++
#include <cstdio>

const int MAXN = 1000005;

long long a[MAXN];
int cnt[63];

int main() {
    int n;
    scanf("%d", &n);
    for (int i = 0; i < n; i++) {
        scanf("%lld", &a[i]);
        for (int j = 0; j < 63; j++) cnt[j] += !!(a[i] & (1ll << j));
    }

    bool isC = true;
    for (int i = 0; i < 63; i++) if (cnt[i] % 3) {
        isC = false;
        break;
    }

    if (isC) {
        puts("C");
    } else {
        bool isA = false;
        for (int i = 0; i < n; i++) {
            bool flag = true;
            long long temp = 0;
            for (int j = 0; j < 63; j++) {
                int c = (cnt[j] - !!(a[i] & (1ll << j))) % 3;
                if (c == 1) {
                    flag = false;
                    break;
                }
                if (c == 2) temp |= (1ll << j);
            }
            if (flag && temp < a[i]) {
                isA = true;
                break;
            }
        }
        puts(isA ? "A" : "B");
    }

    return 0;
}
```
