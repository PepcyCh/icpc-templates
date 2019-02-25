# 树状数组（Binary Indexed Tree/Fenwick Tree）

```c++
#include <cstdio>
#include <algorithm>

const int MAXN = 100005;

struct BIT {
    int n;
    long long a[MAXN];

    static constexpr int lowbit(int x) { return x & -x; }

    void init(int n) {
        this->n = n;
        std::fill(a + 1, a + n + 1, 0);
    }

    void update(int pos, int val) {
        for (int i = pos; i <= n; i += lowbit(i)) a[i] += val;
    }

    long long query(int pos) {
        long long res = 0;
        for (int i = pos; i; i -= lowbit(i)) res += a[i];
        return res;
    }
    long long query(int l, int r) { return query(r) - query(l - 1); }
} bit;

int main() {
    int n;
    scanf("%d", &n);

    static int a[MAXN];
    for (int i = 1; i <= n; i++) scanf("%d", &a[i]);

    bit.init(n);
    for (int i = 1; i <= n; i++) bit.update(i, a[i]);

    int q;
    scanf("%d", &q);
    while (q--) {
        int op;
        scanf("%d", &op);

        if (op == 1) {
            int pos, val;
            scanf("%d %d", &pos, &val);
            bit.update(pos, val);
        } else {
            int l, r;
            scanf("%d %d", &l, &r);
            printf("%lld\n", bit.query(l, r));
        }
    }

    return 0;
}
```
