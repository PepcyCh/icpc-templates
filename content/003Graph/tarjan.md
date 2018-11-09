# Tarjan

三个连通分量与对应的 Tarjan 算法。

### 有向图找强连通分量并缩点建图

```c++
#include <cstdio>
#include <stack>
#include <algorithm>

const int MAXN = 100005;
const int MAXM = 500005;

template <typename T>
struct Edge {
    T *u, *v;
    Edge *next;

    Edge() {}
    Edge(T *u, T *v) : u(u), v(v), next(u->e) {}
};

struct Comp {
    Edge<Comp> *e;
    int size;
} C[MAXN];
Edge<Comp> _poolC[MAXM], *_currC = _poolC;

void addEdge(Comp *u, Comp *v) {
    u->e = new (_currC++) Edge<Comp>(u, v);
}

struct Node {
    Edge<Node> *e;
    Comp *c;
    int dfn, low;
    bool ins;
} N[MAXN];
Edge<Node> _poolN[MAXM], *_currN = _poolN;

void addEdge(int u, int v) {
    N[u].e = new (_currN++) Edge<Node>(&N[u], &N[v]);
}

namespace Tarjan {
    int dfsClock, compCnt;
    std::stack<Node *> s;

    void dfs(Node *u) {
        u->dfn = u->low = ++dfsClock;
        s.push(u);
        u->ins = true;

        for (Edge<Node> *e = u->e; e; e = e->next) {
            if (!e->v->dfn) {
                dfs(e->v);
                u->low = std::min(u->low, e->v->low);
            } else if (e->v->ins) {
                u->low = std::min(u->low, e->v->dfn);
            }
        }

        if (u->dfn == u->low) {
            Comp *c = &C[++compCnt];
            Node *v;

            do {
                v = s.top();
                s.pop();
                v->ins = false;
                v->c = c;
                c->size++;
            } while (u != v);
        }
    }

    void findSCC(int n) {
        dfsClock = compCnt = 0;
        while (!s.empty()) s.pop();

        for (int i = 1; i <= n; i++) if (!N[i].dfn) dfs(&N[i]);
    }
    
    void rebuild() {
        for (Edge<Node> *e = _poolN; e != _currN; e++) if (e->u->c != e->v->c) addEdge(e->u->c, e->v->c);
    }
}

int main() {
    int n, m;
    scanf("%d %d", &n, &m);

    for (int i = 0, u, v; i < m; i++) {
        scanf("%d %d", &u, &v);
        addEdge(u, v);
    }

    Tarjan::findSCC(n);
    Tarjan::rebuild();

    return 0;
}
```

### 无向图找边双连通分量并缩点建图

```c++
#include <cstdio>
#include <vector>
#include <queue>
#include <algorithm>

const int MAXN = 300005;

template <typename T>
struct Edge {
    T *u, *v;
    Edge *next;

    Edge() {}
    Edge(T *u, T *v) : u(u), v(v), next(u->e) {}
};

struct Comp {
    Edge<Comp> *e;
    int dist;
} C[MAXN];
Edge<Comp> _poolC[MAXN << 1], *_currC = _poolC;
void addEdgeC(int u, int v) {
    C[u].e = new (_currC++) Edge<Comp>(&C[u], &C[v]);
    C[v].e = new (_currC++) Edge<Comp>(&C[v], &C[u]);
}

struct Node {
    Edge<Node> *e;
    int dfn, low;
} N[MAXN];
Edge<Node> _poolN[MAXN << 1], *_currN = _poolN;
void addEdgeN(int u, int v) {
    N[u].e = new (_currN++) Edge<Node>(&N[u], &N[v]);
    N[v].e = new (_currN++) Edge<Node>(&N[v], &N[u]);
}

namespace Tarjan {
    struct DJS {
        int f[MAXN];

        void init(int n) {
            for (int i = 1; i <= n; i++) f[i] = i;
        }

        int find(int x) {
            return x == f[x] ? x : f[x] = find(f[x]);
        }

        void merge(int x, int y) {
            f[find(x)] = find(y);
        }
    } djs;

    int dfsClock;
    std::vector<Edge<Node> *> bridges;

    void dfs(Node *u, Node *fa = NULL) {
        u->dfn = ++dfsClock;
        u->low = u->dfn;

        for (Edge<Node> *e = u->e; e; e = e->next) if (e->v != fa) {
            if (e->v->dfn) {
                u->low = std::min(u->low, e->v->dfn);
            } else {
                dfs(e->v, u);
                u->low = std::min(u->low, e->v->low);
                if (e->v->low > u->dfn) {
                    bridges.push_back(e);
                } else {
                    djs.merge(u - N, e->v - N);
                }
            }
        }
    }

    void findECC(int n) {
        djs.init(n);
        dfsClock = 0;
        dfs(&N[1]);
    }

    void rebuild() {
        for (auto e : bridges) {
            int x = djs.find(e->u - N);
            int y = djs.find(e->v - N);
            if (x != y) addEdgeC(x, y);
        }
    }
} // namespace Tarjan

int main() {
    int n, m;
    scanf("%d %d", &n, &m);

    for (int i = 0, u, v; i < m; i++) {
        scanf("%d %d", &u, &v);
        addEdgeN(u, v);
    }

    Tarjan::findECC(n);
    Tarjan::rebuild();

    return 0;
}
```

### 无向图找点双连通分量并缩点建图

```c++
#include <cstdio>
#include <vector>
#include <stack>
#include <algorithm>

const int MAXN = 300005;

template <typename T>
struct Edge {
    T *u, *v;
    Edge *next;

    Edge() {}
    Edge(T *u, T *v) : u(u), v(v), next(u->e) {}
};

struct Comp {
    Edge<Comp> *e;
} C[MAXN];
Edge<Comp> _poolC[MAXN << 1], *_currC = _poolC;

void addEdge(Comp *u, Comp *v) {
    u->e = new (_currC++) Edge<Comp>(u, v);
    v->e = new (_currC++) Edge<Comp>(v, u);
}

struct Node {
    Edge<Node> *e;
    Comp *c; // a cut's comp is meaningless
    int dfn, low;
    bool isCut;
} N[MAXN];
Edge<Node> _poolN[MAXN << 1], *_currN = _poolN;

void addEdge(int u, int v) {
    N[u].e = new (_currN++) Edge<Node>(&N[u], &N[v]);
    N[v].e = new (_currN++) Edge<Node>(&N[v], &N[u]);
}

namespace Tarjan {
    int dfsClock, compCnt;
    std::stack<Edge<Node> *> s;
    std::vector<Node *> cuts;

    void dfs(Node *u, Node *fa = NULL) {
        u->dfn = ++dfsClock;
        u->low = u->dfn;
        int size = 0;

        for (Edge<Node> *e = u->e; e; e = e->next) if (e->v != fa) {
            s.push(e);

            if (e->v->dfn) {
                u->low = std::min(u->low, e->v->dfn);
            } else {
                ++size;
                dfs(e->v, u);
                u->low = std::min(u->low, e->v->low);
                if (e->v->low >= u->dfn) {
                    u->isCut = true;
                    cuts.push_back(u);

                    Comp *c = &C[++compCnt];
                    while (true) {
                        Edge<Node> *t = s.top();
                        s.pop();

                        if (t->u->c != c) t->u->c = c;
                        if (t->v->c != c) t->v->c = c;

                        if (t->u == u && t->v == e->v) break;
                    }
                }
            }
        }

        if (!fa && size > 1) {
            cuts.push_back(u);
            u->isCut = true;
        }
    }

    void findVCC(int n) {
        dfsClock = 0;
        while (!s.empty()) s.pop();
        dfs(&N[1]);
    }

    void rebuild() {
        for (Edge<Node> *e = _poolN; e < _currN; e++) if (e->u->c != e->v->c) addEdge(e->u->c, e->v->c);
    }
}

int main() {
    int n, m;
    scanf("%d %d", &n, &m);

    for (int i = 0, u, v; i < m; i++) {
        scanf("%d %d", &u, &v);
        addEdge(u, v);
    }

    Tarjan::findVCC(n);
    Tarjan::rebuild();

    return 0;
}
```
