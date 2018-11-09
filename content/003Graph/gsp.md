# 全局最短路（Global Shortest Path）

普通的 Floyd、Floyd 传递闭包、Johnson 算法。

### Floyd

```c++
#include <cstdio>
#include <climits>
#include <algorithm>

const int MAXN = 1005;

int dist[MAXN][MAXN];

void floyd(int n) {
    for (int k = 1; k <= n; k++) for (int i = 1; i <= n; i++) if (dist[i][k] < INT_MAX)
        for (int j = 1; j <= n; j++) if (dist[k][j] < INT_MAX)
            dist[i][j] = std::min(dist[i][j], dist[i][k] + dist[k][j]);
}

int main() {
    int n, m;
    scanf("%d %d", &n, &m);

    for (int i = 1; i <= n; i++) std::fill(dist[i] + 1, dist[i] + n + 1, INT_MAX);

    for (int i = 0, u, v, w; i < m; i++) {
        scanf("%d %d %d", &u, &v, &w);
        dist[u][v] = dist[v][u] = std::min(dist[u][v], w);
    }

    floyd(n);
    // dist[i][i] is the length of minimum circle through i

    return 0;
}
```

### Floyd 传递闭包

```c++
#include <cstdio>
#include <bitset>
#include <algorithm>

const int MAXN = 2005;

std::bitset<MAXN> G[MAXN];

void floyd(int n) {
    for (int k = 0; k < n; k++) for (int i = 0; i < n; i++)
        if (G[i][k]) G[i] |= G[k];
}

int main() {
    int n, m;
    scanf("%d %d", &n, &m);

    for (int i = 0, u, v; i < m; i++) {
        scanf("%d %d", &u, &v);
        G[u].set(v);
    }
    for (int i = 0; i < n; i++) G[i].set(i);

    floyd(n);

    return 0;
}
```

### Johnson

用于求稀疏（有向）图的全局最短路。如果边权均非负，可省略第一步的 Bellman-Ford。

```c++
#include <cstdio>
#include <climits>
#include <queue>
#include <algorithm>

const int MAXN = 2005;
const int MAXM = 10005;

struct Edge;
struct Node {
    Edge *e;
    int dist[MAXN], cnt, h; // N[j].dist[i] denotes the shortest path from i to j
    bool inq, vis;
} N[MAXN];

struct Edge {
    Node *u, *v;
    Edge *next;
    int w;

    Edge() {}
    Edge(Node *u, Node *v, int w) : u(u), v(v), w(w), next(u->e) {}
} _pool[MAXM + MAXN], *_curr = _pool;

void addEdge(int u, int v, int w) {
    N[u].e = new (_curr++) Edge(&N[u], &N[v], w);
}

namespace BellmanFord {
    bool bellmanFord(Node *s, int n) {
        std::queue<Node *> q;
        q.push(s);
        s->h = 0;

        while (!q.empty()) {
            Node *u = q.front();
            q.pop();
            u->inq = false;

            for (Edge *e = u->e; e; e = e->next) {
                if (e->v->h > u->h + e->w) {
                    e->v->h = u->h + e->w;

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
    
    void solve(int s, int n) {
        for (int i = 0; i < n; i++) {
            N[i].h = INT_MAX;
            N[i].cnt = 0;
            N[i].inq = false;
        }

        bellmanFord(&N[s], n);
    }
} // namespace BellmanFord

namespace Dijkstra {
    struct HeapNode {
        Node *u;
        int dist;

        HeapNode(int dist, Node *u) : u(u), dist(dist) {}

        bool operator<(const HeapNode &rhs) const {
            return dist > rhs.dist;
        }
    };

    void dijkstra(Node *s, int id) {
        std::priority_queue<HeapNode> q;
        s->dist[id] = 0;
        q.emplace(0, s);

        while (!q.empty()) {
            Node *u = q.top().u;
            q.pop();

            if (u->vis) continue;
            u->vis = true;

            for (Edge *e = u->e; e; e = e->next) {
                if (e->v->dist[id] > u->dist[id] + e->w) {
                    e->v->dist[id] = u->dist[id] + e->w;

                    q.emplace(e->v->dist[id], e->v);
                }
            }
        }
    }

    void solve(int s, int n) {
        for (int i = 1; i <= n; i++) {
            N[i].dist[s] = INT_MAX;
            N[i].vis = false;
        }

        dijkstra(&N[s], s);
    }
} // namespace Dijkstra

int main() {
    int n, m;
    scanf("%d %d", &n, &m);
    for (int i = 0, u, v, w; i < m; i++) {
        scanf("%d %d %d", &u, &v, &w);
        addEdge(u, v, w);
    }

    for (int i = 1; i <= n; i++) addEdge(0, i, 0);
    BellmanFord::solve(0, n + 1);
    for (Edge *e = _pool; e < _curr - n; e++) e->w += e->u->h - e->v->h;
    for (int i = 1; i <= n; i++) {
        Dijkstra::solve(i, n);
        for (int j = 1; j <= n; j++) N[j].dist[i] -= N[i].h - N[j].h;
    }
    
    return 0;
}
```
