# 单源最短路（Single Source Shortest Path）

### Dijkstra

```c++
#include <cstdio>
#include <climits>
#include <queue>
#include <algorithm>

const int MAXN = 100005;
const int MAXM = 500005;

struct Edge;
struct Node;

struct Node {
    Edge *e; // use std::vector<Edge> when dense graph.
    int dist;
    bool vis;
} N[MAXN];

struct Edge {
    Node *u, *v;
    Edge *next;
    int w;

    Edge() {}
    Edge(Node *u, Node *v, int w) : u(u), v(v), w(w), next(u->e) {}
} _pool[MAXM << 1], *_curr = _pool;

void addEdge(int u, int v, int w) {
    N[u].e = new (_curr++) Edge(&N[u], &N[v], w);
    N[v].e = new (_curr++) Edge(&N[v], &N[u], w);
}

namespace Dijkstra {
    struct HeapNode {
        Node *u;
        int dist;

        HeapNode(int dist, Node *u) : u(u), dist(dist) {}

        bool operator<(const HeapNode &rhs) const {
            return dist > rhs.dist;
        }
    };

    void dijkstra(Node *s) {
        std::priority_queue<HeapNode> q;
        s->dist = 0;
        q.emplace(0, s);

        while (!q.empty()) {
            Node *u = q.top().u;
            q.pop();

            if (u->vis) continue;
            u->vis = true;

            for (Edge *e = u->e; e; e = e->next) {
                if (e->v->dist > u->dist + e->w) {
                    e->v->dist = u->dist + e->w;

                    q.emplace(e->v->dist, e->v);
                }
            }
        }
    }

    int solve(int s, int t, int n) {
        for (int i = 1; i <= n; i++) {
            N[i].dist = INT_MAX;
            N[i].vis = false;
        }

        dijkstra(&N[s]);

        return N[t].dist;
    }
}

int main() {
    int n, m, s, t;
    scanf("%d %d %d %d", &n, &m, &s, &t);
    
    for (int i = 0, u, v, w; i < m; i++) {
        scanf("%d %d %d", &u, &v, &w);
        addEdge(u, v, w);
    }

    printf("%d\n", Dijkstra::solve(s, t, n));

    return 0;
}
```

### 分层最短路

可视作最短路 DP。

```c++
#include <cstdio>
#include <climits>
#include <queue>
#include <algorithm>

const int MAXN = 10005;
const int MAXM = 50005;
const int MAXK = 10;

struct Edge;
struct Node;

struct Node {
    Edge *e; // use std::vector<Edeg> when dense graph
    int dist[MAXK];
    bool vis[MAXK];
} N[MAXN];

struct Edge {
    Node *u, *v;
    Edge *next;
    int w;

    Edge() {}
    Edge(Node *u, Node *v, int w) : u(u), v(v), w(w), next(u->e) {}
} _pool[MAXM << 1], *_curr = _pool;

void addEdge(int u, int v, int w) {
    N[u].e = new (_curr++) Edge(&N[u], &N[v], w);
    N[v].e = new (_curr++) Edge(&N[v], &N[u], w);
}

namespace Dijkstra {
    struct HeapNode {
        Node *u;
        int dist, k;

        HeapNode(Node *u, int k, int dist) : u(u), k(k), dist(dist) {}

        bool operator<(const HeapNode &rhs) const {
            return dist > rhs.dist;
        }
    };
    
    void dijkstra(Node *s, int maxK) {
        std::priority_queue<HeapNode> q;
        s->dist[0] = 0;
        q.emplace(s, 0, 0);

        while (!q.empty()) {
            Node *u = q.top().u;
            int k = q.top().k;
            q.pop();

            if (u->vis[k]) continue;
            u->vis[k] = true;

            for (Edge *e = u->e; e; e = e->next) {
                if (e->v->dist[k] > u->dist[k] + e->w) {
                    e->v->dist[k] = u->dist[k] + e->w;

                    q.emplace(e->v, k, e->v->dist[k]);
                }

                if (k < maxK && e->v->dist[k + 1] > u->dist[k]) {
                    e->v->dist[k + 1] = u->dist[k];

                    q.emplace(e->v, k + 1, e->v->dist[k + 1]);
                }
            }
        }
    }

    int solve(int s, int t, int n, int k) {
        for (int i = 1; i <= n; i++) {
            std::fill(N[i].dist, N[i].dist + k + 1, INT_MAX);
            std::fill(N[i].vis, N[i].vis + k + 1, false);
        }

        dijkstra(&N[s], k);

        int res = INT_MAX;
        for (int i = 0; i <= k; i++) res = std::min(res, N[t].dist[i]);
        return res;
    }
}

int main() {
    int n, m, k, s, t;
    scanf("%d %d %d %d %d", &n, &m, &k, &s, &t);

    for (int i = 0, u, v, w; i < m; i++) {
        scanf("%d %d %d", &u, &v, &w);
        addEdge(u, v, w);
    }

    printf("%d\n", Dijkstra::solve(s, t, n, k));
    
    return 0;
}
```

### 队列优化的 Bellman-Ford / SPFA

不应该使用，除非是含负权边的单次最短路（如差分约束）。求负权边的多次最短路时，请使用类似 Johnson 算法的方式。

```c++
#include <cstdio>
#include <climits>
#include <queue>
#include <algorithm>

const int MAXN = 100005;
const int MAXM = 500005;

struct Edge;
struct Node;

struct Node {
    Edge *e; // use std::vetor<Edge> when dense graph
    int dist, cnt;
    bool inq;
} N[MAXN];

struct Edge {
    Node *u, *v;
    Edge *next;
    int w;

    Edge() {}
    Edge(Node *u, Node *v, int w) : u(u), v(v), w(w), next(u->e) {}
} _pool[MAXM << 1], *_curr = _pool;

void addEdge(int u, int v, int w) {
    N[u].e = new (_curr++) Edge(&N[u], &N[v], w);
    N[v].e = new (_curr++) Edge(&N[v], &N[u], w);
}

namespace BellmanFord {
    bool bellmanFord(Node *s, int n) {
        std::queue<Node *> q;
        q.push(s);
        s->dist = 0;

        while (!q.empty()) {
            Node *u = q.front();
            q.pop();
            u->inq = false;

            for (Edge *e = u->e; e; e = e->next) {
                if (e->v->dist > u->dist + e->w) {
                    e->v->dist = u->dist + e->w;

                    if (++e->v->cnt >= n) return false;
                    if (!e->v->inq) {
                        e->v->inq = true;
                        q.push(e->v);
                    }
                }
            }
        }
        return true;
    }
    
    int solve(int s, int t, int n) {
        for (int i = 1; i <= n; i++) {
            N[i].dist = INT_MAX;
            N[i].cnt = 0;
            N[i].inq = false;
        }

        return bellmanFord(&N[s], n) ? N[t].dist : -1;
    }
}

int main() {
    int n, m, s, t;
    scanf("%d %d %d %d", &n, &m, &s, &t);

    for (int i = 0, u, v, w; i < m; i++) {
        scanf("%d %d %d", &u ,&v, &w);
        addEdge(u, v, w);
    }

    printf("%d\n", BellmanFord::solve(s, t, n));

    return 0;
}
```
