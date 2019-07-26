# 稀疏表（Sparse Table）

```c++
#include <cstdio>
#include <algorithm>

const int MAXN = 500005;
const int MAXN_LOG = 20;

struct SparseTable {
    int st[MAXN][MAXN_LOG], log[MAXN];

    void init(int *a, int n) {
        int t = 0;
        for (int i = 0; i <= n; i++) {
            while (1 << (t + 1) <= i) t++;
            log[i] = t;
        }

        for (int i = 1; i <= n; i++) st[i][0] = a[i];

        for (int j = 1; j <= log[n]; j++) {
            for (int i = 1; i <= n; i++) {
                if (i + (1 << (j - 1)) <= n) st[i][j] = std::min(st[i][j - 1], st[i + (1 << (j - 1))][j - 1]);
                else st[i][j] = st[i][j - 1];
            }
        }
    }

    int query(int l, int r) {
        int t = log[r - l + 1];
        return std::min(st[l][t], st[r - (1 << t) + 1][t]);
    }
} st;

int main() {
    int n;
    scanf("%d", &n);

    static int a[MAXN];
    for (int i = 1; i <= n; i++) scanf("%d", &a[i]);

    st.init(a, n);

    int q;
    scanf("%d", &q);
    while (q--) {
        int l, r;
        scanf("%d %d", &l, &r);
        printf("%d\n", st.query(l, r));
    }

    return 0;
}
```
