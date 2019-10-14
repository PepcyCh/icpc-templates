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

ZKW 费用流

```c++
#include <cstdio>
#include <climits>
#include <queue>
#include <algorithm>

const int MAXN = 405;

struct Graph {
    struct Edge {
        int v, cap, flow, cost, rev;

        Edge(int v, int cap, int cost, int rev) : v(v), cap(cap), flow(0), cost(cost), rev(rev) {}
    };

    void addEdge(int u, int v, int cap, int cost) {
        G[u].emplace_back(v, cap, cost, G[v].size());
        G[v].emplace_back(u, 0, -cost, G[u].size() - 1);
    }

    std::vector<Edge> G[MAXN];
    int n;
} G;

class MCMF {
public:
    std::pair<int, int> solve(int s, int t, Graph &G) {
        this->s = s;
        this->t = t;
        this->G = &G;

        int flow = 0, cost = 0;
        int rem = INT_MAX;
        while (argument()) {
            std::fill_n(ptr + 1, G.n, 0);
            int d = dfs(s, rem - flow);
            flow += d;
            cost += d * dist[t];
        }

        return std::make_pair(flow, cost);
    }

private:
    int s, t;
    Graph *G;
    int dist[MAXN], ptr[MAXN];
    bool vis[MAXN];

    bool argument() {
        std::fill_n(dist + 1, G->n, INT_MAX);
        std::fill_n(vis + 1, G->n, false);
        dist[s] = 0;

        std::queue<int> q;
        q.push(s);
        while (!q.empty()) {
            int u = q.front();
            q.pop();

            vis[u] = false;
            for (auto &e : G->G[u]) {
                if (e.cap > e.flow && dist[e.v] > dist[u] + e.cost) {
                    dist[e.v] = dist[u] + e.cost;
                    if (!vis[e.v]) {
                        vis[e.v] = true;
                        q.push(e.v);
                    }
                }
            }
        }

        return dist[t] != INT_MAX;
    }

    int dfs(int u, int r) {
        if (u == t) return r;

        vis[u] = true;
        int res = 0;

        for (int &i = ptr[u]; i < G->G[u].size(); i++) {
            auto &e = G->G[u][i];
            if (e.cap > e.flow && dist[e.v] == dist[u] + e.cost && !vis[e.v]) {
                int d = dfs(e.v, std::min(r - res, e.cap - e.flow));
                res += d;
                e.flow += d;
                G->G[e.v][e.rev].flow -= d;
                if (res == r) {
                    vis[u] = false;
                    break;
                }
            }
        }

        return res;
    }
} mcmf;

int main() {
    int n, m;
    scanf("%d %d", &n, &m);

    G.n = n;
    for (int i = 0, u, v, c, w; i < m; i++) {
        scanf("%d %d %d %d", &u, &v, &c, &w);
        G.addEdge(u, v, c, w);
    }

    auto [flow, cost] = mcmf.solve(1, n, G);
    printf("%d %d\n", flow, cost);

    return 0;
}
```
