# 最小树形图

朱-刘算法，O(nm)。

```c++
#include <cstdio>
#include <climits>
#include <algorithm>

const int MAXN = 3005;

namespace ChuLiu {
    int G[MAXN][MAXN], n, eg[MAXN], queue[MAXN];
    bool used[MAXN], pass[MAXN], more;

    void combine(int id, long long &res) {
        int tot = 0, from;
        for (; id != 0 && !pass[id]; id = eg[id]) {
            queue[tot++] = id;
            pass[id] = true;
        }
        for (from = 0; from < tot && queue[from] != id; ++from) {}
        if (from == tot) return;
        more = true;
        for (int i = from; i < tot; i++) {
            res += G[eg[queue[i]]][queue[i]];
            if (i != from) {
                used[queue[i]] = true;
                for (int j = 1; j <= n; j++) if (!used[j])
                    G[id][j] = std::min(G[id][j], G[queue[i]][j]);
            }
        }
        for (int i = 1; i <= n; i++) if (!used[i] && i != id) {
            for (int j = from; j < tot; j++) {
                int k = queue[j];
                G[i][id] = std::min(G[i][id], G[i][k] - G[eg[k]][k]);
            }
        }
    }

    long long solve(int n, int root) {
        ChuLiu::n = n;
        std::fill(used + 1, used + n + 1, false);
        long long res = 0;
        for (more = true; more; ) {
            more = false;
            std::fill(eg + 1, eg + n + 1, 0);
            for (int i = 1, k; i <= n; i++) if (!used[i] && i != root) {
                k = 0;
                for (int j = 1; j <= n; j++) if (!used[j] && i != j) {
                    if (k == 0 || G[j][i] < G[k][i])
                        k = j;
                }
                eg[i] = k;
            }
            std::fill(pass + 1, pass + n + 1, false);
            for (int i = 1; i <= n; i++) if (!used[i] && !pass[i] && i != root)
                combine(i, res);
        }
        for (int i = 1; i <= n; i++) if (!used[i] && i != root) {
            if (G[eg[i]][i] == INT_MAX) return -1;
            res += G[eg[i]][i];
        }
        return res;
    }
}

int main() {
    int n, m, rt;
    scanf("%d %d %d", &n, &m, &rt);
    for (int i = 1; i <= n; i++)
        std::fill(ChuLiu::G[i] + 1, ChuLiu::G[i] + n + 1, INT_MAX);
    for (int i = 0, u, v, w; i < m; i++) {
        scanf("%d %d %d", &u, &v, &w);
        ChuLiu::G[u][v] = w;
    }

    long long ans = ChuLiu::solve(n, rt);
    printf("%lld\n", ans);

    return 0;
}
```

Tarjan 的优化实现，O(n^2)。

```c++
#include <cstdio>
#include <vector>
#include <stack>
#include <queue>
#include <algorithm>

/*
 * It takes O(n^2) time on dense graph
 * To get O(m \log n) time on sparse graph, change 'std::vector<Edge *> in' to __gnu_pbds::priority_queue<Edge *>,
 * and change line 71-72, 113-138
 */

const int MAXN = 2505;

struct Graph {
    struct Edge {
        int u, v, w, ow;
        bool removed;
        std::vector<Edge *> ch;
        Edge *pa;

        Edge() {}
        Edge(int u, int v, int w) : u(u), v(v), w(w), ow(w), pa(NULL), removed(false) {}
    };
    std::vector<Edge> E;

    void addEdge(int u, int v, int w) {
        if (u == v) return;
        E[eid] = Edge(u, v, w);
        G[u].push_back(&E[eid]);
        in[v].push_back(&E[eid]);
        ++eid;
    }

    std::vector<Edge *> G[MAXN], in[MAXN];
    int n, eid;

    void init(int n, int m) {
        this->n = n;
        E.resize(m);
        this->eid = 0;
    }
} G;

struct DJS {
    int f[MAXN];

    void init(int n) {
        for (int i = 0; i <= n; i++) f[i] = i;
    }
    int find(int x) { return x == f[x] ? x : f[x] = find(f[x]); }
    bool test(int x, int y) { return find(x) == find(y); }
    void merge(int x, int y) { f[find(y)] = find(x); }
};

class OptimumBranching {
public:
    long long solve(Graph &G, int rt) {
        std::fill_n(change + 1, G.n, 0);
        S.init(G.n);
        W.init(G.n);

        std::stack<int> roots;
        std::vector<Edge *> F;

        for (int i = 1; i <= G.n; i++) if (i != rt) roots.push(i);
        while (!roots.empty()) {
            int u = roots.top();
            roots.pop();

            Edge *min = NULL;
            for (Edge *e : G.in[u]) if (!min || e->w < min->w) min = e;

            F.push_back(min);
            for (Edge *e : cycle[u]) {
                e->pa = min;
                min->ch.push_back(e);
            }

            if (cycle[u].empty()) lambda[u] = min;

            if (!W.test(min->u, min->v)) {
                enter[u] = min;
                W.merge(min->u, min->v);
            } else {
                std::vector<Edge *> cycleEdge;
                std::vector<int> cycleRepr;

                Edge *max = min;
                enter[u] = NULL;

                cycleEdge.push_back(min);
                cycleRepr.push_back(S.find(min->v));
                for (int v = S.find(min->u); enter[v]; v = S.find(enter[v]->u)) {
                    cycleEdge.push_back(enter[v]);
                    cycleRepr.push_back(v);
                    if (max->w < enter[v]->w) max = enter[v];
                }

                for (Edge *e : cycleEdge) change[S.find(e->v)] = max->w - e->w;

                int nr = cycleRepr.front();
                for (int v : cycleRepr) {
                    S.merge(v, nr);
                    nr = S.find(nr);
                }
                roots.push(nr);
                cycle[nr].swap(cycleEdge);

                for (int v : cycleRepr) for (Edge *e : G.in[v]) e->w += change[v];

                std::vector<Edge *> nin;
                for (int i = 1; i < cycleRepr.size(); i++) {
                    auto i1 = G.in[cycleRepr[i]].begin();
                    auto e1 = G.in[cycleRepr[i]].end();
                    auto i2 = G.in[cycleRepr[i - 1]].begin();
                    auto e2 = G.in[cycleRepr[i - 1]].end();

                    while (i1 != e1 || i2 != e2) {
                        while (i1 != e1 && S.test((*i1)->u, nr)) ++i1;
                        while (i2 != e2 && S.test((*i2)->u, nr)) ++i2;

                        if (i1 == e1 && i2 == e2) break;
                        else if (i1 == e1) nin.push_back(*i2++);
                        else if (i2 == e2) nin.push_back(*i1++);
                        else if ((*i1)->u < (*i2)->u) nin.push_back(*i1++);
                        else if ((*i1)->u > (*i2)->u) nin.push_back(*i2++);
                        else {
                            if ((*i1)->w < (*i2)->w) nin.push_back(*i1);
                            else nin.push_back(*i2);
                            ++i1;
                            ++i2;
                        }
                    }

                    G.in[cycleRepr[i]].swap(nin);
                    nin.clear();
                }
                G.in[nr].swap(G.in[cycleRepr.back()]);
                change[nr] = 0;
            }
        }

        long long res = 0;
        std::stack<Edge *> froots;
        for (Edge *e : F) if (e->pa == NULL) froots.push(e);
        while (!froots.empty()) {
            Edge *e = froots.top();
            froots.pop();
            if (e->removed) continue;
            res += e->ow;
            remove(lambda[e->v], froots);
        }

        return res;
    }

private:
    using Edge = Graph::Edge;
    DJS S, W;
    std::vector<Edge *> cycle[MAXN];
    Edge *lambda[MAXN], *enter[MAXN];
    int change[MAXN];

    void remove(Edge *e, std::stack<Edge *> &roots) {
        for (; e; e = e->pa) {
            e->removed = true;
            for (Edge *c : e->ch) {
                roots.push(c);
                c->pa = NULL;
            }
            e->ch.clear();
        }
    }
} ob;

bool check(int rt, int n) {
    static int vis[MAXN];

    std::queue<int> q;
    q.push(rt);
    vis[rt] = true;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (Graph::Edge *e : G.G[u]) if (!vis[e->v]) {
            vis[e->v] = true;
            q.push(e->v);
        }
    }
    for (int i = 1; i <= n; i++) if (!vis[i]) return false;
    return true;
}

int main() {
    int n, m, rt;
    scanf("%d %d %d", &n, &m, &rt);

    G.init(n, m);
    for (int i = 0, u, v, w; i < m; i++) {
        scanf("%d %d %d", &u, &v, &w);
        G.addEdge(u, v, w);
    }

    if (!check(rt, n)) return puts("-1"), 0;

    long long ans = ob.solve(G, rt);
    printf("%lld\n", ans);

    return 0;
}
```
