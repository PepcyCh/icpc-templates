# 并查集（Disjoint Set）

```c++
#include <cstdio>

const int MAXN = 1000005;

struct DisjointSet {
    int fa[MAXN];

    void init(int n) {
        for (int i = 0; i < n; i++) fa[i] = i;
    }

    int find(int x) {
        return x == fa[x] ? x : fa[x] = find(fa[x]);
    }

    void merge(int x, int y) {
        fa[find(y)] = find(x);
    }

    bool test(int x, int y) {
        return find(x) == find(y);
    }
} djs;

int main() {
    int n, q;
    scanf("%d %d", &n, &q);

    djs.init(n);

    while (q--) {
        int op, u, v;
        scanf("%d %d %d", &op, &u, &v);
        if (op == 0) djs.merge(u, v);
        else puts(djs.test(u, v) ? "Yes" : "No");
    }

    return 0;
}
```
