# 网络流

Dinic 与 ISAP，以及最小割输出方案。

### Dinic

```c++
#include <cstdio>
#include <climits>
#include <queue>
#include <vector>
#include <algorithm>

const int MAXN = 1000005;
const int MAXM = 4000005;

struct Edge;
struct Node;

struct Node {
    std::vector<Edge> e; // use Edge* when sparse graph
    Edge *curr;
    int level;
} N[MAXN];

struct Edge {
    Node *u,*v;
    int cap, flow, rev;

    Edge(Node *u, Node *v, int cap, int rev) : u(u), v(v), cap(cap), flow(0), rev(rev) {}
};

void addEdge(int u, int v, int cap) {
    N[u].e.emplace_back(&N[u], &N[v], cap, N[v].e.size());
    N[v].e.emplace_back(&N[v], &N[u], 0, N[u].e.size() - 1);
}

namespace Dinic {
    bool level(Node *s, Node *t, int n) {
        for (int i = 0; i < n; i++) N[i].level = 0;
        
        static std::queue<Node *> q;
        while (!q.empty()) q.pop();
        q.push(s);
        s->level = 1;

        while (!q.empty()) {
            Node *u = q.front();
            q.pop();

            for (Edge *e = &u->e.front(); e <= &u->e.back(); e++) {
                if (e->cap > e->flow && e->v->level == 0) {
                    e->v->level = u->level + 1;
                    if (e->v == t) return true;
                    q.push(e->v);
                }
            }
        }

        return false;
    }

    int findPath(Node *u, Node *t, int limit = INT_MAX) {
        if (u == t) return limit;

        for (Edge *&e = u->curr; e <= &u->e.back(); e++) {
            if (e->cap > e->flow && e->v->level == u->level + 1) {
                int flow = findPath(e->v, t, std::min(limit, e->cap - e->flow));
                if (flow > 0) {
                    e->flow += flow;
                    e->v->e[e->rev].flow -= flow;
                    return flow;
                }
            }
        }

        return 0;
    }

    int solve(int s, int t, int n) {
        int res = 0;

        while (level(&N[s], &N[t], n)) {
            for (int i = 0; i < n; i++) N[i].curr = &N[i].e.front();
            int flow;
            while ((flow = findPath(&N[s], &N[t])) > 0) res += flow;
        }

        return res;
    }
}

int main() {
    int n, m, s, t;
    scanf("%d %d %d %d", &n, &m, &s, &t);

    for (int i = 0, u, v, w; i < m; i++) {
        scanf("%d %d %d", &u, &v, &w);
        addEdge(u, v, w);
    }

    printf("%d\n", Dinic::solve(s, t, n));

    return 0;
}
```

### ISAP

```c++
#include <cstdio>
#include <climits>
#include <algorithm>

const int MAXN = 105;
const int MAXM = 5005;

struct Edge;
struct Node {
    Edge *e, *curr; // use std::vector<Edge> when dense graph
    int dist;
} N[MAXN];

struct Edge {
    Node *u, *v;
    Edge *next, *rev;
    int cap, flow;

    Edge() {}
    Edge(Node *u, Node *v, int cap) : u(u), v(v), cap(cap), flow(0), next(u->e) {}
} _pool[MAXM << 1], *_curr = _pool;

void addEdge(int u, int v, int cap) {
    N[u].e = new (_curr++) Edge(&N[u], &N[v], cap);
    N[v].e = new (_curr++) Edge(&N[v], &N[u], 0);
    N[u].e->rev = N[v].e;
    N[v].e->rev = N[u].e;
}

namespace ISAP {
    int cnt[MAXN], n;
    Node *s, *t;

    int flow(Node *u, int limit = INT_MAX) {
        if (u == t) return limit;

        int temp = limit;
        for (Edge *&e = u->curr; e; e = e->next) if (u->dist == e->v->dist + 1 && e->cap > e->flow) {
            int f = flow(e->v, std::min(temp, e->cap - e->flow));
            e->flow += f;
            e->rev->flow -= f;
            temp -= f;
            if (!temp) return limit;
        }

        if (!(--cnt[u->dist++])) s->dist = n + 1;
        ++cnt[u->dist];
        u->curr = u->e;

        return limit - temp;
    }

    long long solve(int s, int t, int n) {
        ISAP::n = n;
        ISAP::s = &N[s];
        ISAP::t = &N[t];
        for (int i = 1; i <= n; i++) {
            N[i].curr = N[i].e;
            N[i].dist = 0;
            cnt[i] = 0;
        }
        cnt[0] = n;

        long long res = 0;
        while (N[s].dist <= n) res += flow(&N[s]);
        return res;
    }
}

int main() {
    int n, m, s, t;
    scanf("%d %d %d %d", &n, &m, &s, &t);

    for (int i = 0, u, v, c; i < m; i++) {
        scanf("%d %d %d", &u, &v, &c);
        addEdge(u, v, c);
    }

    long long ans = ISAP::solve(s, t, n);
    printf("%lld\n", ans);

    return 0;
}
```

### 最小割输出方案

```c++
std::vector<Edge> cuts;
void bfs(int s, int t, int n) {
    static std::queue<Node *> q;
    while (!q.empty()) q.pop();
    q.push(&N[s]);
    N[s].vis = true;

    while (!q.empty()) {
        Node *u = q.front();
        q.pop();

        for (Edge e : u->e) if (e.cap > e.flow && !e.v->vis) {
            e.v->vis = true;
            q.push(e.v);
        }
    }

    for (int i = 0; i < n; i++) if (N[i].vis) for (Edge e : N[i].e) {
        if (e.cap == e.flow && e.cap > 0 && !e.v->vis) {
            cust.push_back(e);
        }
    }
}
```
