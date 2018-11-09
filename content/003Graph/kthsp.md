# 两点间 k 短路

```c++
#include <cstdio>
#include <climits>
#include <queue>
#include <algorithm>

const int MAXN = 1005;
const int MAXM = 100005;

struct Edge;
struct Node {
    Edge *e, *eo;
    int dist, times;
    bool vis;
} N[MAXN];

struct Edge {
    Node *u, *v;
    Edge *next;
    int w;

    Edge() {}
    Edge(Node *u, Node *v, int w, Edge *next) : u(u), v(v), w(w), next(next) {}
} _pool[MAXM << 1], *_curr = _pool;

void addEdge(int u, int v, int w) {
    N[u].e = new (_curr++) Edge(&N[u], &N[v], w, N[u].e);
    N[v].eo = new (_curr++) Edge(&N[v], &N[u], w, N[v].eo);
}

namespace Dijkstra {
    struct HeapNode {
        Node *u;
        int dist;

        HeapNode(Node *u, int dist) : u(u), dist(dist) {}

        bool operator<(const HeapNode &rhs) const {
            return dist > rhs.dist;
        }
    };

    void dijkstra(Node *s) {
        static std::priority_queue<HeapNode> q;
        s->dist = 0;
        q.emplace(s, 0);
        while (!q.empty()) {
            Node *u = q.top().u;
            q.pop();

            if (u->vis) continue;
            u->vis = true;

            for (Edge *e = u->eo; e; e = e->next) {
                if (e->v->dist > u->dist + e->w) {
                    e->v->dist = u->dist + e->w;
                    q.emplace(e->v, e->v->dist);
                }
            }
        }
    }

    void solve(int s, int n) {
        for (int i = 1; i <= n; i++) {
            N[i].dist = INT_MAX;
            N[i].vis = false;
        }

        dijkstra(&N[s]);
    }
}

namespace KthShortest {
    struct HeapNode {
        Node *u;
        int curr, last;

        HeapNode(Node *u, int curr, int last) : u(u), curr(curr), last(last) {}

        bool operator<(const HeapNode &rhs) const {
            return curr + last > rhs.curr + rhs.last;
        }
    };

    int astar(Node *s, Node *t, int k) {
        static std::priority_queue<HeapNode> q;
        q.emplace(s, 0, s->dist);
        while (!q.empty()) {
            Node *u = q.top().u;
            int curr = q.top().curr;
            int last = q.top().last;
            q.pop();

            ++u->times;
            if (u->times == k && u == t) return curr + last;
            if (u->times > k) continue;

            for (Edge *e = u->e; e; e = e->next) q.emplace(e->v, curr + e->w, e->v->dist);
        }

        return -1;
    }

    int solve(int s, int t, int k, int n) {
        Dijkstra::solve(t, n);
        return astar(&N[s], &N[t], k);
    }
}

int main() {
    int n, m;
    scanf("%d %d", &n, &m);
    for (int i = 0, u, v, w; i < m; i++) {
        scanf("%d %d %d", &u, &v, &w);
        addEdge(u, v, w);
    }

    int s, t, k;
    scanf("%d %d %d", &s, &t, &k);
    if (s == t) ++k;

    int ans = KthShortest::solve(s, t, k, n);
    printf("%d\n", ans);

    return 0;
}
```
