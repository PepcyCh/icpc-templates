# 一般图最大匹配

带花树（Blossom Algorithm）。

```c++
#include <bits/stdc++.h>

const int MAXN = 505;

struct Graph {
    std::vector<int> G[MAXN];
    int n, mate[MAXN];

    void addEdge(int u, int v) {
        G[u].push_back(v);
        G[v].push_back(u);
    }
} G;

class Blossom {
  public:
    int solve(Graph &G) {
        this->G = &G;
        for (int i = 1; i <= G.n; i++) G.mate[i] = -1;
        for (int i = 1; i <= G.n; i++) if (G.mate[i] == -1) aug(i);

        int res = 0;
        for (int i = 1; i <= G.n; i++) res += (G.mate[i] > i);
        return res;
    }

  private:
    struct DJS {
        int f[MAXN];

        void init(int n) {
            for (int i = 1; i <= n; i++) f[i] = i;
        }

        int find(int x) { return x == f[x] ? x : f[x] = find(f[x]); }
        bool test(int x, int y) { return find(x) == find(y); }
        void merge(int x, int y) { f[find(x)] = find(y); }
    } djs;

    int next[MAXN], dsu[MAXN], mark[MAXN], vis[MAXN];
    Graph *G;

    int lca(int x, int y) {
        static int t = 0;
        ++t;
        for (; ; std::swap(x, y)) if (x != -1) {
            if (vis[x = djs.find(x)] == t) return x;
            vis[x] = t;
            x = (G->mate[x] != -1) ? next[G->mate[x]] : -1;
        }
    }

    std::queue<int> q;
    void group(int a, int p) {
        for (int b, c; a != p; djs.merge(a, b), djs.merge(b, c), a = c) {
            b = G->mate[a], c = next[b];
            if (djs.find(c) != p) next[c] = b;
            if (mark[b] == 2) mark[b] = 1, q.push(b);
            if (mark[c] == 2) mark[c] = 1, q.push(c);
        }
    }

    void aug(int s) {
        for (int i = 1; i <= G->n; i++) {
            next[i] = vis[i] = -1;
            mark[i] = 0;
        }
        djs.init(G->n);
        while (!q.empty()) q.pop();

        q.push(s);
        mark[s] = 1;
        while (G->mate[s] == -1 && !q.empty()) {
            int x = q.front();
            q.pop();

            for (int y : G->G[x]) {
                if (y != G->mate[x] && !djs.test(x, y) && mark[y] != 2) {
                    if (mark[y] == 1) {
                        int p = lca(x, y);
                        if (djs.find(x) != p) next[x] = y;
                        if (djs.find(y) != p) next[y] = x;
                        group(x, p);
                        group(y, p);
                    } else if (G->mate[y] == -1) {
                        next[y] = x;
                        for (int j = y, k, l; j != -1; j = l) {
                            k = next[j];
                            l = G->mate[k];
                            G->mate[j] = k;
                            G->mate[k] = j;
                        }
                        break;
                    } else {
                        next[y] = x;
                        q.push(G->mate[y]);
                        mark[G->mate[y]] = 1;
                        mark[y] = 2;
                    }
                }
            }
        }
    }
} blossom;

int main() {
    int n, m;
    scanf("%d %d", &n, &m);

    G.n = n;
    for (int i = 0, u, v; i < m; i++) {
        scanf("%d %d", &u, &v);
        G.addEdge(u, v);
    }

    blossom.solve(G);

    for (int i = 1; i <= n; i++) printf("%d%c", G.mate[i], " \n"[i == n]);
    
    return 0;
}
```
