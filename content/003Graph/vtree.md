# 虚树

```c++
#include <cstdio>
#include <algorithm>

const int MAXN = 300005;
const int MAXN_LOG = 20;

template <typename T>
struct Edge {
    T *u, *v;
    Edge *next;

    Edge() {}
    Edge(T *u, T *v) : u(u), v(v), next(u->e) {}
};

struct Node {
    Edge<Node> *e;
    Node *f[MAXN_LOG];
    int dfn, dep;
} N[MAXN];
Edge<Node> _poolN[MAXN << 1], *_currN = _poolN;
void addEdge(int u, int v) {
    N[u].e = new (_currN++) Edge<Node>(&N[u], &N[v]);
    N[v].e = new (_currN++) Edge<Node>(&N[v], &N[u]);
}

struct VirN {
    Edge<VirN> *e;
    Node *r;
} V[MAXN];
Edge<VirN> _poolV[MAXN << 1], *_currV = _poolV;
void addEdge(VirN *u, VirN *v) {
    u->e = new (_currV++) Edge<VirN>(u, v);
    v->e = new (_currV++) Edge<VirN>(v, u);
}

void dfs(Node *u, bool init = true) {
    if (init) {
        u->dep = 1;
        u->f[0] = u;
    }

    static int dfsClock = 0;
    u->dfn = ++dfsClock;

    for (int i = 1; i < MAXN_LOG; i++) u->f[i] = u->f[i - 1]->f[i - 1];

    for (Edge<Node> *e = u->e; e; e = e->next) if (e->v != u->f[0]) {
        e->v->f[0] = u;
        e->v->dep = u->dep + 1;
        dfs(e->v, false);
    }
}

Node *lca(Node *u, Node *v) {
    if (u->dep < v->dep) std::swap(u, v);

    for (int i = MAXN_LOG - 1; ~i; i--) {
        if (u->f[i]->dep >= v->dep) u = u->f[i];
    }

    for (int i = MAXN_LOG - 1; ~i; i--) {
        if (u->f[i] != v->f[i]) {
            u = u->f[i];
            v = v->f[i];
        }
    }

    return u == v ? u : u->f[0];
}

void build(bool flag, int n, int &tot) {
    static VirN *stack[MAXN];
    int top = 0;

    if (!flag) stack[top++] = &V[0];

    for (int i = 1; i <= n; i++) {
        if (!top) {
            stack[top++] = &V[i];
            continue;
        }
        Node *p = lca(stack[top - 1]->r, V[i].r);
        if (p == stack[top - 1]->r) {
            stack[top++] = &V[i];
            continue;
        }
        while (top - 2 >= 0 && stack[top - 2]->r->dep >= p->dep) {
            addEdge(stack[top - 2], stack[top - 1]);
            top--;
        }
        if (stack[top - 1]->r != p) {
            V[++tot].r = p;
            addEdge(&V[tot], stack[--top]);
            stack[top++] = &V[tot];
        }
        stack[top++] = &V[i];
    }

    for (int i = 0; i < top - 1; i++) addEdge(stack[i], stack[i + 1]);
}

bool cmp(int i, int j) { return N[i].dfn < N[j].dfn; }

void clear(int n) {
    _currV = _poolV;
    for (int i = 0; i <= n; i++) V[i].e = NULL;
}

void solve() {
    int k;
    scanf("%d", &k);
    clear(k);

    static int order[MAXN], h[MAXN];
    for (int i = 0; i < k; i++) scanf("%d", &h[i]), order[i] = h[i];

    std::sort(h, h + k, cmp);
    int tot = 0;
    bool flag = h[0] == 1;
    for (int i = 0; i < k; i++) V[++tot].r = &N[h[i]];
    build(flag, tot, tot);

    // do something...
}

int main() {
    int n;
    scanf("%d", &n);
    for (int i = 1, u, v; i < n; i++) {
        scanf("%d %d", &u, &v);
        addEdge(u, v);
    }

    dfs(&N[1]);

    V[0].r = &N[1];
    int q;
    scanf("%d", &q);
    while (q--) solve();

    return 0;
}
```
