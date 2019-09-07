# 2-SAT

```c++
#include <cstdio>
#include <stack>
#include <algorithm>

const int MAXN = 1005;
const int MAXM = 2000005;

template<typename T>
struct Edge {
    T *u, *v;
    Edge *next;

    Edge() {}
    Edge(T *u, T *v) : u(u), v(v), next(u->e) {}
};

struct Comp {
    Edge<Comp> *e;
    Comp *opp;
    int deg, mark;

    Comp() : mark(-1) {}
} C[MAXN << 1];
Edge<Comp> _poolC[MAXM], *_currC = _poolC;
void addEdge(Comp *u, Comp *v) {
    u->e = new (_currC++) Edge<Comp>(u, v);
    ++v->deg;
}

struct Node {
    Edge<Node> *e;
    int dfn, low;
    Comp *belong;
    bool ins;
} N[MAXN << 1];
Edge<Node> _poolN[MAXM], *_currN = _poolN;
void addEdge(int u, int v) {
    N[u].e = new (_currN++) Edge<Node>(&N[u], &N[v]);
}

void rebuild() {
    for (Edge<Node> *e = _poolN; e != _currN; e++) if (e->u->belong != e->v->belong)
        addEdge(e->v->belong, e->u->belong);
}

int getV(int x, int k) {
    return (x << 1) - k;
}

namespace TwoSat {
    std::stack<Node *> s;
    int dfsClock, sccCnt;

    void tarjan(Node *u) {
        s.push(u);
        u->dfn = u->low = ++dfsClock;
        u->ins = true;

        for (Edge<Node> *e = u->e; e; e = e->next) {
            if (e->v->dfn == 0) {
                tarjan(e->v);
                u->low = std::min(u->low, e->v->low);
            } else if (e->v->ins) {
                u->low = std::min(u->low, e->v->dfn);
            }
        }

        if (u->dfn == u->low) {
            Node *v;
            Comp *c = &C[sccCnt++];
            do {
                v = s.top();
                s.pop();
                v->ins = false;
                v->belong = c;
            } while (v != u);
        }
    }

    bool check(int n) {
        for (int i = 1; i <= n; i++) if (N[getV(i, 0)].belong == N[getV(i, 1)].belong)
            return false;
        return true;
    }

    bool solve(int n) {
        sccCnt = dfsClock = 0;
        for (int i = 1; i <= n << 1; i++) if (N[i].dfn == 0) tarjan(&N[i]);
        return check(n);
    }
} // namespace TwoSat

namespace TopoSort {
    std::stack<Comp *> s;

    void dfs(Comp *u) {
        if (u->mark != -1) return;

        u->mark = 0;
        for (Edge<Comp> *e; e; e = e->next) dfs(e->v);
    }

    void solve(int n) {
        for (int i = 1; i <= n; i++) if (C[i].deg == 0) s.push(&C[i]);

        while (!s.empty()) {
            Comp *u = s.top();
            s.pop();

            if (u->mark != -1) continue;

            u->mark = 1;
            dfs(u->opp);

            for (Edge<Comp> *e = u->e; e; e = e->next) if (!(--e->v->deg))
                s.push(e->v);
        }
    }
} // namespace TopoSort

int main() {
    int n, m;
    scanf("%d %d", &n, &m);

    for (int i = 0, u, v, op; i < m; i++) {
        scanf("%d %d %d", &u, &v, &op);

        if (op == 1) { // u and v can't be true at the same time
            addEdge(getV(u, 1), getV(v, 0));
            addEdge(getV(v, 1), getV(u, 0));
        } else if (op == 2) { // u and v can't be false at the same time
            addEdge(getV(u, 0), getV(v, 1));
            addEdge(getV(v, 0), getV(u, 1));
        } else if (op == 3) { // u and v must be different
            addEdge(getV(u, 1), getV(v, 0));
            addEdge(getV(v, 1), getV(u, 0));
            addEdge(getV(u, 0), getV(v, 1));
            addEdge(getV(v, 0), getV(u, 1));
        } else { // u and v must be the same
            addEdge(getV(u, 1), getV(v, 1));
            addEdge(getV(v, 1), getV(u, 1));
            addEdge(getV(u, 0), getV(v, 0));
            addEdge(getV(v, 0), getV(u, 0));
        }
    }

    TwoSat::solve(n);

    rebuild();

    TopoSort::solve(TwoSat::sccCnt);

    return 0;
}
```
