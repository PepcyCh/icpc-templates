# 斯坦纳树（Steiner Tree）

```c++
#include <cstdio>
#include <climits>
#include <queue>
#include <algorithm>

const int MAXN = 50005;
const int MAXM = 300005;
const int MAXD = 5;

struct Edge;
struct Node {
    Edge *e;
    bool vis[1 << MAXD];
    long long f[1 << MAXD];
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
        int s;
        long long dist;

        HeapNode(Node *u, int s, long long dist) : u(u), s(s), dist(dist) {}

        bool operator<(const HeapNode &rhs) const {
            return dist > rhs.dist;
        }
    };

    std::priority_queue<HeapNode> q;
    void dijkstra() {
        while (!q.empty()) {
            Node *u = q.top().u;
            int s = q.top().s;
            q.pop();

            if (u->vis[s]) continue;
            u->vis[s] = true;

            for (Edge *e = u->e; e; e = e->next) {
                if (e->v->f[s] > u->f[s] + e->w) {
                    e->v->f[s] = u->f[s] + e->w;
                    q.emplace(e->v, s, e->v->f[s]);
                }
            }
        }
    }

    void init(int n, int d) {
        for (int i = 1; i <= n; i++) for (int s = 0; s < 1 << d; s++) {
            N[i].f[s] = LLONG_MAX >> 1ll;
            N[i].vis[s] = false;
        }
        while (!q.empty()) q.pop();
    }
}

// O(3^n + 2^n * m \log n)
long long steiner(int n, int *p, int d) {
    Dijkstra::init(n, d);
    for (int i = 0; i < d; i++) N[p[i]].f[1 << i] = 0;

    for (int S = 0; S < 1 << d; S++) {
        for (int s = (S - 1) & S; s; s = (s - 1) & S)
            for (int i = 1; i <= n; i++) N[i].f[S] = std::min(N[i].f[S], N[i].f[s] + N[i].f[S ^ s]);
        for (int i = 1; i <= n; i++) if (N[i].f[S] < LLONG_MAX >> 1ll) Dijkstra::q.emplace(&N[i], S, N[i].f[S]);
        Dijkstra::dijkstra();
    }

    long long res = LLONG_MAX;
    for (int i = 1; i <= n; i++) res = std::min(res, N[i].f[(1 << d) - 1]);
    return res;
}

int main() {
    int n, m, d;
    scanf("%d %d %d", &n, &m, &d);
    for (int i = 0, u, v, w; i < m; i++) {
        scanf("%d %d %d", &u, &v, &w);
        addEdge(u, v, w);
    }

    static int p[MAXN];
    for (int i = 0; i < d; i++) scanf("%d", &p[i]);

    long long ans = steiner(n, p, d);
    printf("%lld\n", ans);

    return 0;
}
```
