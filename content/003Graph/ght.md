# Gomory-Hu Tree / 最小割树

```c++
#include <cstdio>
#include <climits>
#include <vector>
#include <queue>
#include <algorithm>

const int MAXN = 3005;
const int MAXN_LOG = 13;

struct Edge;
struct Node {
    std::vector<Edge> e;
    Edge *curr;
    int level;
    bool vis;
} N[MAXN];

struct Edge {
    Node *u, *v;
    int cap, flow, rev;

    Edge(Node *u, Node *v, int cap, int rev) : u(u), v(v), cap(cap), flow(0), rev(rev) {}
};

bool G[MAXN][MAXN];

void addEdge(int u, int v, int cap) {
    N[u].e.emplace_back(&N[u], &N[v], cap, N[v].e.size());
    N[v].e.emplace_back(&N[v], &N[u], cap, N[u].e.size() - 1);
    G[u][v] = G[v][u] = true;
}

namespace Dinic {
    bool level(Node *s, Node *t, int n) {
        for (int i = 1; i <= n; i++) N[i].level = 0;
        static std::queue<Node *> q;
        q.push(s);
        s->level = 1;
        while (!q.empty()) {
            Node *u = q.front();
            q.pop();

            for (Edge &e : u->e) {
                if (e.cap > e.flow && e.v->level == 0) {
                    e.v->level = u->level + 1;
                    q.push(e.v);
                }
            }
        }
        return t->level;
    }

    int findPath(Node *u, Node *t, int limit = INT_MAX) {
        if (u == t) return limit;
        int res = 0;
        for (Edge *&e = u->curr; e && e <= &u->e.back(); e++) {
            if (e->cap > e->flow && e->v->level == u->level + 1) {
                int flow = findPath(e->v, t, std::min(limit, e->cap - e->flow));
                if (flow > 0) {
                    e->flow += flow;
                    e->v->e[e->rev].flow -= flow;
                    limit -= flow;
                    res += flow;
                    if (!limit) return res;
                } else e->v->level = -1;
            }
        }
        return res;
    }

    int solve(int s, int t, int n) {
        for (int i = 1; i <= n; i++) for (Edge &e : N[i].e) e.flow = 0;

        int res = 0;
        while (level(&N[s], &N[t], n)) {
            for (int i = 1; i <= n; i++) N[i].curr = &N[i].e.front();
            int flow;
            while ((flow = findPath(&N[s], &N[t])) > 0) res += flow;
        }
        return res;
    }

    void bfs(Node *s, int n) {
        for (int i = 1; i <= n; i++) N[i].vis = false;

        static std::queue<Node *> q;
        q.push(s);
        s->vis = true;
        while (!q.empty()) {
            Node *u = q.front();
            q.pop();
            for (const Edge &e : u->e) {
                if (e.cap > e.flow && !e.v->vis) {
                    e.v->vis = true;
                    q.push(e.v);
                }
            }
        }
    }
}

namespace GomoryHuTree {
    struct Edge;
    struct Node {
        Edge *e;
        Node *f[MAXN_LOG];
        int dep, min[MAXN_LOG];
    } N[MAXN];

    struct Edge {
        Node *u, *v;
        Edge *next;
        int w;

        Edge() {}
        Edge(Node *u, Node *v, int w) : u(u), v(v), next(u->e), w(w) {}
    } _pool[MAXN << 1], *_curr;

    void addEdge(int u, int v, int w) {
        N[u].e = new (_curr++) Edge(&N[u], &N[v], w);
        N[v].e = new (_curr++) Edge(&N[v], &N[u], w);
    }

    void dfs(Node *u, Node *fa = NULL) {
        u->f[0] = (fa ? fa : u);
        u->dep = (fa ? fa->dep : 0) + 1;
        for (int i = 1; i < MAXN_LOG; i++) {
            u->f[i] = u->f[i - 1]->f[i - 1];
            u->min[i] = std::min(u->min[i - 1], u->f[i - 1]->min[i - 1]);
        }

        for (Edge *e = u->e; e; e = e->next) if (e->v != fa) {
            e->v->min[0] = e->w;
            dfs(e->v, u);
        }
    }

    int query(Node *u, Node *v) {
        if (u->dep < v->dep) std::swap(u, v);
        int res = INT_MAX;

        for (int i = MAXN_LOG - 1; ~i; i--) if (u->f[i]->dep >= v->dep) {
            res = std::min(res, u->min[i]);
            u = u->f[i];
        }

        for (int i = MAXN_LOG - 1; ~i; i--) if (u->f[i] != v->f[i]) {
            res = std::min({res, u->min[i], v->min[i]});
            u = u->f[i];
            v = v->f[i];
        }

        if (u != v) res = std::min({res, u->min[0], v->min[0]});

        return res;
    }
    int query(int u, int v) { return query(&N[u], &N[v]); }

    void build(const std::vector<int> &nodes, int n) {
        if (nodes.size() <= 1) return;
        int s = nodes[0], t = nodes[1];
        int flow = Dinic::solve(s, t, n);
        addEdge(s, t, flow);

        Dinic::bfs(&::N[s], n);

        std::vector<int> ln, rn;
        for (int u : nodes) {
            if (::N[u].vis) ln.push_back(u);
            else rn.push_back(u);
        }

        build(ln, n);
        build(rn, n);
    }

    void build(int n) {
        _curr = _pool;
        std::vector<int> vec(n);
        for (int i = 1; i <= n; i++) vec[i - 1] = i;
        build(vec, n);

        N[1].min[0] = INT_MAX;
        dfs(&N[1]);
    }
}

int main() {
    int n, m;
    scanf("%d %d", &n, &m);
    for (int i = 0, u, v, w; i < m; i++) {
        scanf("%d %d %d", &u, &v, &w);
        addEdge(u, v, w);
    }

    GomoryHuTree::build(n);

    int q;
    scanf("%d", &q);
    while (q--) {
        int u, v;
        scanf("%d %d", &u, &v);

        int ans = GomoryHuTree::query(u, v);
        printf("%d\n", ans);
    }

    return 0;
}
```
