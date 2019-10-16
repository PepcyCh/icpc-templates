# DSU on Tree / 树上启发式合并

```c++
#include <cstdio>
#include <vector>
#include <algorithm>

const int MAXN = 100005;

struct Graph {
    void addEdge(int u, int v) {
        G[u].push_back(v);
        G[v].push_back(u);
    }

    std::vector<int> G[MAXN];
    long long ans[MAXN];
    int col[MAXN];
    int n;
} G;

class DSUT {
public:
    void solve(Graph &G) {
        this->G = &G;
        findHeavy(1);
        dfs(1);
    }

private:
    Graph *G;
    int cnt[MAXN];
    int size[MAXN], son[MAXN];
    int skip, max;
    long long sum;

    void findHeavy(int u, int fa = -1) {
        size[u] = 1;
        for (int v : G->G[u]) if (v != fa) {
            findHeavy(v, u);
            size[u] += size[v];
            if (size[v] > size[son[u]]) son[u] = v;
        }
    }

    void edit(int u, int fa, int d) {
        cnt[G->col[u]] += d;
        if (d > 0) {
            if (cnt[G->col[u]] > max) {
                sum = G->col[u];
                max = cnt[G->col[u]];
            } else if (cnt[G->col[u]] == max) {
                sum += G->col[u];
            }
        }
        for (int v : G->G[u]) if (v != fa && v != skip) edit(v, u, d);
    }

    void dfs(int u, int fa = -1, bool keep = false) {
        for (int v : G->G[u]) if (v != fa && v != son[u]) dfs(v, u);
        if (son[u]) {
            dfs(son[u], u, true);
            skip = son[u];
        }
        edit(u, fa, 1);
        G->ans[u] = sum;
        skip = -1;
        if (!keep) {
            edit(u, fa, -1);
            sum = 0;
            max = 0;
        }
    }
} dsut;

int main() {
    int n;
    scanf("%d", &n);
    G.n = n;

    for (int i = 1; i <= n; i++) scanf("%d", &G.col[i]);
    for (int i = 1, u, v; i < n; i++) {
        scanf("%d %d", &u, &v);
        G.addEdge(u, v);
    }

    dsut.solve(G);

    for (int i = 1; i <= n; i++) printf("%lld%c", G.ans[i], " \n"[i == n]);

    return 0;
}
```
