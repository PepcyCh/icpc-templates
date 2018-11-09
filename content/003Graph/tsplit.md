# 树链剖分

### 重链剖分

```c++
#include <cstdio>
#include <climits>
#include <algorithm>

const int MAXN = 100005;

struct Node;
struct Edge;
struct Chain;

struct Node {
    Edge *e;
    Chain *c;
    Node *max, *fa;
    int dfn, dfnR, dep, size, val; // use dfnR when there're operations on subtrees
} N[MAXN];

struct Edge {
    Node *u, *v;
    Edge *next;

    Edge() {}
    Edge(Node *u, Node *v) : u(u), v(v), next(u->e) {}
} _poolE[MAXN << 1], *_currE = _poolE;
void addEdge(int u, int v) {
    N[u].e = new (_currE++) Edge(&N[u], &N[v]);
    N[v].e = new (_currE++) Edge(&N[v], &N[u]);
}

struct Chain {
    Node *top;

    Chain(Node *top = NULL) : top(top) {}
} _poolC[MAXN], *_currC = _poolC;

void dfs1(Node *u) {
    u->size = 1;

    for (Edge *e = u->e; e; e = e->next) if (e->v != u->fa) {
        e->v->dep = u->dep + 1;
        e->v->fa = u;

        dfs1(e->v);

        u->size += e->v->size;
        if (!u->max || u->max->size < e->v->size) u->max = e->v;
    }
}
void dfs2(Node *u, Node *fa = NULL) {
    static int dfsClock = 0;
    u->dfn = ++dfsClock;

    if (!u->fa || u->fa->max != u) u->c = new (_currC++) Chain(u);
    else u->c = u->fa->c;

    if (u->max) dfs2(u->max);
    for (Edge *e = u->e; e; e = e->next) if (e->v != u->fa && e->v != u->max) dfs2(e->v);

    u->dfnR = dfsClock;
}
void split() {
    N[1].dep = 1;
    dfs1(&N[1]);
    dfs2(&N[1]);
}

struct SegT {
    struct Node {
        int l, r;
        Node *lc, *rc;
        int max, tag;

        Node() {}
        Node(int pos, int val) : l(pos), r(pos), max(val), tag(0), lc(NULL), rc(NULL) {}
        Node(Node *lc, Node *rc) : l(lc->l), r(rc->r), lc(lc), rc(rc), tag(0) {
            maintain();
        }

        void add(int d) {
            max += d;
            tag += d;
        }

        void pushDown() {
            if (tag) {
                lc->add(tag);
                rc->add(tag);
                tag = 0;
            }
        }

        void maintain() {
            max = std::max(lc->max, rc->max);
        }

        void update(int l, int r, int d) {
            if (l > this->r || this->l > r) return;
            if (l <= this->l && this->r <= r) {
                add(d);
                return;
            }
            pushDown();
            lc->update(l, r, d);
            rc->update(l, r, d);
            maintain();
        }

        int query(int l, int r) {
            if (l > this->r || this->l > r) return INT_MIN;
            if (l <= this->l && this->r <= r) return max;
            pushDown();
            return std::max(lc->query(l, r), rc->query(l, r));
        }
    } *root, _pool[MAXN << 1], *_curr;

    SegT() : root(NULL), _curr(_pool) {}
    
    Node *_build(int l, int r, int *a) {
        if (l == r) return new (_curr++) Node(l, a[l]);
        int mid = l + ((r - l) >> 1);
        return new (_curr++) Node(_build(l, mid, a), _build(mid + 1, r, a));
    }
    void build(int l, int r, int *a) {
        root = _build(l, r, a);
    }

    void update(int l, int r, int d) {
        root->update(l, r, d);
    }

    int query(int l, int r) {
        return root->query(l, r);
    }
} segT;

void update(int a, int b, int d) {
    Node *u = &N[a], *v = &N[b];

    while (u->c != v->c) {
        if (u->c->top->dep < v->c->top->dep) std::swap(u, v);
        segT.update(u->c->top->dfn, u->dfn, d);
        u = u->c->top->fa;
    }

    if (u->dep > v->dep) std::swap(u, v);
    segT.update(u->dfn, v->dfn, d);
}

int query(int a, int b) {
    Node *u = &N[a], *v = &N[b];
    int res = INT_MIN;

    while (u->c != v->c) {
        if (u->c->top->dep < v->c->top->dep) std::swap(u, v);
        res = std::max(res, segT.query(u->c->top->dfn, u->dfn));
        u = u->c->top->fa;
    }

    if (u->dep > v->dep) std::swap(u, v);
    res = std::max(res, segT.query(u->dfn, v->dfn));

    return res;
}

int main() {
    int n;
    scanf("%d", &n);

    for (int i = 1; i <= n; i++) scanf("%d", &N[i].val);

    for (int i = 1, u, v; i < n; i++) {
        scanf("%d %d", &u, &v);
        addEdge(u, v);
    }

    split();

    static int temp[MAXN];
    for (int i = 1; i <= n; i++) temp[N[i].dfn] = N[i].val;
    segT.build(1, n, temp);

    int q;
    scanf("%d", &q);

    while (q--) {
        int op, u, v;
        scanf("%d %d %d", &op, &u, &v);
        if (op == 1) {
            int d;
            scanf("%d", &d);
            update(u, v, d);
        } else printf("%d\n", query(u, v));
    }

    return 0;
}
```
