# 上下界网络流

无源无汇上下界可行流、有源有汇上细节最大流、有源有汇上下界最小流；其中有源有汇上下界最大流可改造为上下界费用流。

### 无源无汇上下界可行流
```c++
#include <cstdio>
#include <climits>
#include <vector>
#include <queue>
#include <algorithm>

const int MAXN = 205;
const int MAXM = 10205;

struct Edge;
struct Node;

struct Node {
    std::vector<Edge> e;
    Edge *curr;
    int level, extra;
} N[MAXN];

struct Edge {
    Node *u, *v;
    int cap, flow, rev;

    Edge(Node *u, Node *v, int cap, int rev)
        : u(u), v(v), cap(cap), flow(0), rev(rev) {}
};

struct Pair {
    int u, num, lower;

    Pair() {}
    Pair(int u, int num, int lower) : u(u), num(num), lower(lower) {}

    int flow() const {
        return N[u].e[num].flow + lower;
    }
} E[MAXM];

Pair addEdge(int u, int v, int lower, int upper) {
    int cap = upper - lower;
    N[u].e.push_back(Edge(&N[u], &N[v], cap, N[v].e.size()));
    N[v].e.push_back(Edge(&N[v], &N[u], 0, N[u].e.size() - 1));

    N[u].extra -= lower;
    N[v].extra += lower;

    return Pair(u, N[u].e.size() - 1, lower);
}

namespace Dinic {
    bool level(Node *s, Node *t, int n) {
        for (int i = 0; i < n; i++) N[i].level = 0;
        static std::queue<Node *> q;
        q.push(s);
        s->level = 1;
        while (!q.empty()) {
            Node *u = q.front();
            q.pop();
            for (Edge *e = &u->e.front(); e <= &u->e.back(); e++) {
                if (e->cap > e->flow && e->v->level == 0) {
                    e->v->level = u->level + 1;
                    q.push(e->v);
                }
            }
        }
        return t->level;
    }

    int findPath(Node *u, Node *t, int limit = INT_MAX) {
        if (u == t) return limit;
        int res = 0;
        for (Edge *&e = u->curr; e <= &u->e.back(); e++) {
            if (e->cap > e->flow && e->v->level == u->level + 1) {
                int flow = findPath(e->v, t, std::min(limit, e->cap - e->flow));
                if (flow > 0) {
                    e->flow += flow;
                    e->v->e[e->rev].flow -= flow;
                    limit -= flow;
                    res += flow;
                    if (limit <= 0) return res;
                } else e->v->level = -1;
            }
        }
        return res;
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
    int n, m;
    scanf("%d %d", &n, &m);
    const int s = 0, t = n + 1;

    for (int i = 0, u, v, lower, upper; i < m; i++) {
        scanf("%d %d %d %d", &u, &v, &lower, &upper);
        E[i] = addEdge(u, v, lower, upper);
    }

    int sum = 0;
    for (int i = 1; i <= n; i++) {
        if (N[i].extra > 0) {
            sum += N[i].extra;
            addEdge(s, i, 0, N[i].extra);
        } else if (N[i].extra < 0) {
            addEdge(i, t, 0, -N[i].extra);
        }
    }

    int maxFlow = Dinic::solve(s, t, n + 2);
    if (maxFlow < sum) {
        puts("NO");
    } else {
        puts("YES");
        for (int i = 0; i < m; i++) printf("%d\n", E[i].flow());
    }

    return 0;
}
```

### 有源有汇上下界最大流

若要求「有源有汇上下界费用流」，将网络流改为 Edmonds-Karp 即可，最后费用是两次的费用和。

```c++
#include <cstdio>
#include <climits>
#include <vector>
#include <queue>
#include <algorithm>

const int MAXN = 205;
const int MAXM = 10000;

struct Edge;
struct Node;

struct Node {
    std::vector<Edge> e;
    Edge *curr;
    int level, extra;
} N[MAXN];

struct Edge {
    Node *u, *v;
    int cap, flow, rev;

    Edge(Node *u, Node *v, int cap, int rev) : u(u), v(v), cap(cap), flow(0), rev(rev) {}
};

void addEdge(int u, int v, int lower, int upper) {
    int cap = upper - lower;
    N[u].e.push_back(Edge(&N[u], &N[v], cap, N[v].e.size()));
    N[v].e.push_back(Edge(&N[v], &N[u], 0, N[u].e.size() - 1));

    N[u].extra -= lower;
    N[v].extra += lower;
}

namespace Dinic {
    bool level(Node *s, Node *t, int n) {
        for (int i = 0; i < n; i++) N[i].level = 0;
        static std::queue<Node *> q;
        q.push(s);
        s->level = 1;
        while (!q.empty()) {
            Node *u = q.front();
            q.pop();
            for (Edge *e = &u->e.front(); e <= &u->e.back(); e++) {
                if (e->cap > e->flow && e->v->level == 0) {
                    e->v->level = u->level + 1;
                    q.push(e->v);
                }
            }
        }
        return t->level;
    }

    int findPath(Node *u, Node *t, int limit = INT_MAX) {
        if (u == t) return limit;
        int res = 0;
        for (Edge *&e = u->curr; e <= &u->e.back(); e++) {
            if (e->cap > e->flow && e->v->level == u->level + 1) {
                int flow = findPath(e->v, t, std::min(limit, e->cap - e->flow));
                if (flow > 0) {
                    e->flow += flow;
                    e->v->e[e->rev].flow -= flow;
                    limit -= flow;
                    res += flow;
                    if (limit <= 0) return res;
                } else e->v->level = -1;
            }
        }
        return res;
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
    const int S = 0, T = n + 1;

    for (int i = 0, u, v, lower, upper; i < m; i++) {
        scanf("%d %d %d %d", &u, &v, &lower, &upper);
        addEdge(u, v, lower, upper);
    }
    addEdge(t, s, 0, INT_MAX);
    int sum = 0;
    for (int i = 1; i <= n; i++) {
        if (N[i].extra > 0) {
            sum += N[i].extra;
            addEdge(S, i, 0, N[i].extra);
        } else if (N[i].extra < 0) {
            addEdge(i, T, 0, -N[i].extra);
        }
    }

    int maxFlow = Dinic::solve(S, T, n + 2);
    if (maxFlow < sum) {
        puts("No");
    } else {
        puts("Yes");
        printf("%d\n", Dinic::solve(s, t, n + 2));
    }

    return 0;
}
```

### 有源有汇上下界最小流

```c++
#include <cstdio>
#include <climits>
#include <vector>
#include <queue>
#include <algorithm>

const int MAXN = 50005;
const int MAXM = 125005;

struct Edge;
struct Node;

struct Node {
    Edge *e, *curr;
    int level, extra;
} N[MAXN];

struct Edge {
    Node *u, *v;
    Edge *next, *rev;
    int cap, flow;

    Edge() {}
    Edge(Node *u, Node *v, int cap) : u(u), v(v), cap(cap), flow(0), next(u->e) {}
} _pool[MAXM + MAXN << 1], *_curr = _pool;

Edge *addEdge(int u, int v, int lower, int upper) {
    int cap = upper - lower;
    N[u].e = new (_curr++) Edge(&N[u], &N[v], cap);
    N[v].e = new (_curr++) Edge(&N[v], &N[u], 0);
    (N[u].e->rev = N[v].e)->rev = N[u].e;

    N[u].extra -= lower;
    N[v].extra += lower;

    return N[u].e;
}

namespace Dinic {
    bool level(Node *s, Node *t, int n) {
        for (int i = 0; i < n; i++) N[i].level = 0;
        static std::queue<Node *> q;
        q.push(s);
        s->level = 1;
        while (!q.empty()) {
            Node *u = q.front();
            q.pop();
            for (Edge *e = u->e; e; e = e->next) {
                if (e->cap > e->flow && e->v->level == 0) {
                    e->v->level = u->level + 1;
                    q.push(e->v);
                }
            }
        }
        return t->level;
    }

    int findPath(Node *u, Node *t, int limit = INT_MAX) {
        if (u == t) return limit;
        int res = 0;
        for (Edge *&e = u->curr; e; e = e->next) {
            if (e->cap > e->flow && e->v->level == u->level + 1) {
                int flow = findPath(e->v, t, std::min(limit, e->cap - e->flow));
                if (flow > 0) {
                    e->flow += flow;
                    e->rev->flow -= flow;
                    limit -= flow;
                    res += flow;
                    if (limit <= 0) return res;
                } else e->v->level = -1;
            }
        }
        return res;
    }

    int solve(int s, int t, int n) {
        int res = 0;
        while (level(&N[s], &N[t], n)) {
            for (int i = 0; i < n; i++) N[i].curr = N[i].e;
            int flow;
            while ((flow = findPath(&N[s], &N[t])) > 0) res += flow;
        }
        return res;
    }
}

int main() {
    int n, m, s, t;
    scanf("%d %d %d %d", &n, &m, &s, &t);
    const int S = 0, T = n + 1;

    for (int i = 0, u, v, lower, upper; i < m; i++) {
        scanf("%d %d %d %d", &u, &v, &lower, &upper);
        addEdge(u, v, lower, upper);
    }
    Edge *e = addEdge(t, s, 0, INT_MAX);
    int sum = 0;
    for (int i = 1; i <= n; i++) {
        if (N[i].extra > 0) {
            sum += N[i].extra;
            addEdge(S, i, 0, N[i].extra);
        } else if (N[i].extra < 0) {
            addEdge(i, T, 0, -N[i].extra);
        }
    }

    int maxFlow = Dinic::solve(S, T, n + 2);
    if (maxFlow < sum) {
        puts("No");
    } else {
        puts("Yes");
        int flow = e->flow;
        e->cap = e->rev->cap = 0;
        printf("%d\n", flow - Dinic::solve(t, s, n + 2));
    }

    return 0;
}
```
