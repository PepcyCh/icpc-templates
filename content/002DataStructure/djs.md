# 并查集

```c++
#include <cstdio>

const int MAXN = 1000005;

struct DSU {
    int fa[MAXN];

    void init(int n) { for (int i = 0; i < n; i++) fa[i] = i; }
    int find(int x) { return x == fa[x] ? x : fa[x] = find(fa[x]); }
    void merge(int x, int y) { fa[find(y)] = find(x); }
    bool test(int x, int y) { return find(x) == find(y); }
} dsu;

int main() {
    int n, q;
    scanf("%d %d", &n, &q);

    dsu.init(n);

    while (q--) {
        int op, u, v;
        scanf("%d %d %d", &op, &u, &v);
        if (op == 0) dsu.merge(u, v);
        else puts(dsu.test(u, v) ? "Yes" : "No");
    }

    return 0;
}
```
