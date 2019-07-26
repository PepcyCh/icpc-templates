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

const int MAXN = 2505;
const int MAXM = MAXN * MAXN;

struct Edge;
struct Node {
    std::vector<Edge *> e, in, cycle; // if no need to traverse, e is not needed
    Edge *lambda, *enter;
    int change;
    bool vis;
} N[MAXN];

struct Edge {
    Node *u, *v;
    std::vector<Edge *> ch;
    Edge *pa;
    int w, ow;
    bool removed;

    Edge() {}
    Edge(Node *u, Node *v, int w) : u(u), v(v), w(w), ow(w), pa(NULL), removed(false) {}
} E[MAXM];
int _curr;

void addEdge(int u, int v, int w) {
    if (u == v) return;
    E[_curr] = Edge(&N[u], &N[v], w);
    N[v].in.push_back(&E[_curr]);
    N[u].e.push_back(&E[_curr]);
    _curr++;
}

namespace OptimumBranching {
    struct DJS {
        int f[MAXN];

        void init(int n) {
            for (int i = 0; i <= n; i++) f[i] = i;
        }
        int find(int x) { return x == f[x] ? x : f[x] = find(f[x]); }
        bool test(int a, int b) { return find(a) == find(b); }
        void merge(int a, int b) { f[find(b)] = find(a); }
    } S, W;

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

    long long solve(int n, int rt) {
        for (int i = 1; i <= n; i++) N[i].change = 0;
        S.init(n), W.init(n);

        static std::stack<Node *> roots;
        static std::vector<Edge *> F;

        for (int i = 1; i <= n; i++) if (i != rt) roots.push(&N[i]);

        while (!roots.empty()) {
            Node *u = roots.top();
            roots.pop();

            Edge *min = NULL;
            for (Edge *e : u->in) if (min == NULL || e->w < min->w) min = e;

            F.push_back(min);
            for (Edge *e : u->cycle) {
                e->pa = min;
                min->ch.push_back(e);
            }

            if (u->cycle.empty()) u->lambda = min;

            if (!W.test(min->u - N, min->v - N)) {
                u->enter = min;
                W.merge(min->u - N, min->v - N);
            } else {
                static std::vector<Edge *> cycle;
                static std::vector<Node *> cycleRepr;
                cycle.clear(), cycleRepr.clear();
                Edge *max = min;
                u->enter = NULL;

                cycle.push_back(min);
                cycleRepr.push_back(&N[S.find(min->v - N)]);
                for (Node *v = &N[S.find(min->u - N)]; v->enter; v = &N[S.find(v->enter->u - N)]) {
                    cycle.push_back(v->enter);
                    cycleRepr.push_back(v);
                    if (max->w < v->enter->w) max = v->enter;
                }

                for (Edge *e : cycle) N[S.find(e->v - N)].change = max->w - e->w;

                Node *nr = cycleRepr.front();
                for (Node *v : cycleRepr) {
                    S.merge(v - N, nr - N);
                    nr = &N[S.find(nr - N)];
                }
                roots.push(nr);
                nr->cycle.swap(cycle);

                for (Node *v : cycleRepr) for (Edge *e : v->in) e->w += v->change;

                static std::vector<Edge *> nin;
                for (int i = 1; i < cycleRepr.size(); i++) {
                    auto i1 = cycleRepr[i]->in.begin();
                    auto e1 = cycleRepr[i]->in.end();
                    auto i2 = cycleRepr[i - 1]->in.begin();
                    auto e2 = cycleRepr[i - 1]->in.end();

                    while (i1 != e1 || i2 != e2) {
                        while (i1 != e1 && S.test((*i1)->u - N, nr - N)) ++i1;
                        while (i2 != e2 && S.test((*i2)->u - N, nr - N)) ++i2;

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

                    cycleRepr[i]->in.swap(nin);
                    nin.clear();
                }
                nr->in.swap(cycleRepr.back()->in);
                nr->change = 0;
            }
        }

        long long res = 0;
        static std::stack<Edge *> froots;
        for (Edge *e : F) if (e->pa == NULL) froots.push(e);
        while (!froots.empty()) {
            Edge *e = froots.top();
            froots.pop();
            if (e->removed) continue;
            res += e->ow;
            remove(e->v->lambda, froots);
        }

        return res;
    }
}

bool check(int rt, int n) {
    static std::queue<Node *> q;
    q.push(&N[rt]);
    N[rt].vis = true;
    while (!q.empty()) {
        Node *u = q.front();
        q.pop();
        for (Edge *e : u->e) if (!e->v->vis) {
            e->v->vis = true;
            q.push(e->v);
        }
    }
    for (int i = 1; i <= n; i++) if (!N[i].vis) return false;
    return true;
}

int main() {
    int n, m, rt;
    scanf("%d %d %d", &n, &m, &rt);
    for (int i = 0, u, v, w; i < m; i++) {
        scanf("%d %d %d", &u, &v, &w);
        addEdge(u, v, w);
    }

    if (!check(rt, n)) return puts("-1"), 0;

    long long ans = OptimumBranching::solve(n, rt);
    printf("%lld\n", ans);

    return 0;
}
```
