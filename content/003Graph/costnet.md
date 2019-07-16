# 费用流

Dijkstra 费用流

```c++
#include <cstdio>
#include <climits>
#include <vector>
#include <queue>
#include <algorithm>

const int MAXN = 405;

struct Edge;
struct Node;

struct Node {
    std::vector<Edge> e;
    Edge *pre;
    int flow, h, dist;
} N[MAXN];

struct Edge {
    Node *u, *v;
    int cap, flow, cost, rev;

    Edge(Node *u, Node *v, int cap, int cost, int rev) : u(u), v(v), rev(rev), cap(cap), flow(0), cost(cost) {}
};

void addEdge(int u, int v, int cap, int cost) {
    N[u].e.push_back(Edge(&N[u], &N[v], cap, cost, N[v].e.size()));
    N[v].e.push_back(Edge(&N[v], &N[u], 0, -cost, N[u].e.size() - 1));
}

namespace EdmondsKarp {
    struct HeapNode {
        Node *u;
        int dist;

        HeapNode(Node *u, int dist) : u(u), dist(dist) {}

        bool operator<(const HeapNode &rhs) const {
            return dist > rhs.dist;
        }
    };

    void solve(int s, int t, int n, int &flow, int &cost) {
        flow = cost = 0;
        // if there exists negative cost, run bellman-ford on h[] first
        while (true) {
            for (int i = 1; i <= n; i++) {
                N[i].dist = INT_MAX;
                N[i].flow = 0;
                N[i].pre = NULL;
            }

            std::priority_queue<HeapNode> q;
            q.push(HeapNode(&N[s], 0));

            N[s].dist = 0;
            N[s].flow = INT_MAX;

            while (!q.empty()) {
                HeapNode un = q.top();
                q.pop();
                Node *u = un.u;

                if (u->dist != un.dist) continue;

                for (Edge *e = &u->e.front(); e <= &u->e.back(); e++) {
                    int newCost = e->cost + u->h - e->v->h;
                    if (e->cap > e->flow && e->v->dist > u->dist + newCost) {
                        e->v->dist = u->dist + newCost;
                        e->v->flow = std::min(u->flow, e->cap - e->flow);
                        e->v->pre = e;

                        q.push(HeapNode(e->v, e->v->dist));
                    }
                }
            }

            if (N[t].dist == INT_MAX) break; // minimum cost maximum flow
            // if (N[t].dist + N[t].h > 0) break; // minimum cost available flow

            for (int i = 1; i <= n; i++) N[i].h = std::min(N[i].h + N[i].dist, INT_MAX >> 1);

            for (Edge *e = N[t].pre; e; e = e->u->pre) {
                e->flow += N[t].flow;
                e->v->e[e->rev].flow -= N[t].flow;
            }

            flow += N[t].flow;
            cost += N[t].h * N[t].flow;
        }
    }
}

int main() {
    int n, m;
    scanf("%d %d", &n, &m);

    for (int i = 0, u, v, c, w; i < m; i++) {
        scanf("%d %d %d %d", &u, &v, &c, &w);
        addEdge(u, v, c, w);
    }

    int flow, cost;
    EdmondsKarp::solve(1, n, n, flow, cost);
    printf("%d %d\n", flow, cost);

    return 0;
}
```
