# 主席树/可持久化线段树（Persistent Segment Tree）

主席树求区间 k 小值和树上两点间 k 小值。

### 区间 k 小值

```c++
#include <cstdio>
#include <algorithm>

const int MAXN = 200005;
const int MAXN_LOG = 20;

template <typename T, size_t SIZE>
struct MemoryPool {
    char mem[sizeof(T) * SIZE], *top;

    MemoryPool() : top(mem) {}

    void *alloc() {
        char *res = top;
        top += sizeof (T);
        return (void *) res;
    }
};

class PSegT {
  private:
    struct Node {
        Node *lc, *rc;
        int cnt;
        static MemoryPool<Node, MAXN * MAXN_LOG> pool;
        static Node *nil;

        Node(int cnt) : lc(nil), rc(nil), cnt(cnt) {
            static bool init = true;
            if (init) { lc = rc = this; init = false; }
        }
        Node(Node *lc = nil, Node *rc = nil) : lc(lc), rc(rc), cnt(lc->cnt + rc->cnt) {}

        void *operator new(size_t) {
            return pool.alloc();
        }

        Node *insert(int l, int r, int val) {
            if (val == l && val == r) return new Node(cnt + 1);
            int mid = l + ((r - l) >> 1);
            if (val <= mid) return new Node(lc->insert(l, mid, val), rc);
            else return new Node(lc, rc->insert(mid + 1, r, val));
        }

        int rank() const { return lc->cnt; }
    } *root[MAXN];
    int n;

  public:
    void build(int *a, int n) {
        this->n = n;
        root[0] = new Node();
        for (int i = 1; i <= n; i++) root[i] = root[i - 1]->insert(1, n, a[i]);
    }

    int query(int l, int r, int k) {
        Node *L = root[l - 1], *R = root[r];
        int min = 1, max = n;
        while (min < max) {
            int mid = min + ((max - min) >> 1), temp = R->rank() - L->rank();
            if (k <= temp) L = L->lc, R = R->lc, max = mid;
            else L = L->rc, R = R->rc, k -= temp, min = mid + 1;
        }
        return min;
    }
} pSegT;
MemoryPool<PSegT::Node, MAXN * MAXN_LOG> PSegT::Node::pool;
PSegT::Node *PSegT::Node::nil = new PSegT::Node(0);

int map[MAXN];
void discrete(int *a, int n) {
    std::copy(a + 1, a + n + 1, map);
    std::sort(map, map + n);
    int *end = std::unique(map, map + n);
    for (int i = 1; i <= n; i++)
        a[i] = std::lower_bound(map, end, a[i]) - map + 1;
}

int main() {
    int n, q;
    scanf("%d %d", &n, &q);

    static int a[MAXN];
    for (int i = 1; i <= n; i++) scanf("%d", &a[i]);
    discrete(a, n);

    pSegT.build(a, n);

    while (q--) {
        int l, r, k;
        scanf("%d %d %d", &l, &r, &k);
        printf("%d\n", map[pSegT.query(l, r, k) - 1]);
    }

    return 0;
}
```

### 树上两点间 k 小值

```c++
#include <cstdio>
#include <climits>
#include <algorithm>

const int MAXN = 100005;
const int MAXN_LOG = 18;

template <typename T, size_t SIZE>
struct MemoryPool {
    char mem[sizeof(T) * SIZE], *top;

    MemoryPool() : top(mem) {}

    void *alloc() {
        char *res = top;
        top += sizeof (T);
        return (void *) res;
    }
};

struct PSegT *null;
struct PSegT {
    PSegT *lc, *rc;
    int cnt;
    static MemoryPool<PSegT, MAXN * 40> pool;

    PSegT(int cnt) : cnt(cnt), lc(null), rc(null) {}
    PSegT(PSegT *lc, PSegT *rc) : lc(lc), rc(rc), cnt(lc->cnt + rc->cnt) {}

    void *operator new(size_t) {
        return pool.alloc();
    }

    PSegT *insert(int l, int r, int val) {
        if (val < l || val > r) return null;
        if (val == l && val == r) return new PSegT(cnt + 1);
        int mid = l + (r - l) / 2;
        if (val <= mid) return new PSegT(lc->insert(l, mid, val), rc);
        else return new PSegT(lc, rc->insert(mid + 1, r, val));
    }
};
MemoryPool<PSegT, MAXN * 40> PSegT::pool;

struct Edge;
struct Node;

struct Node {
    Edge *e;
    PSegT *seg;
    Node *f[MAXN_LOG];
    int dep, val;
} N[MAXN];

struct Edge {
    Node *u, *v;
    Edge *next;

    Edge() {}
    Edge(Node *u, Node *v) : u(u), v(v), next(u->e) {}
} _pool[MAXN << 1], *_curr = _pool;

void addEdge(int u, int v) {
    N[u].e = new (_curr++) Edge(&N[u], &N[v]);
    N[v].e = new (_curr++) Edge(&N[v], &N[u]);
}

void init() {
    null = new PSegT(0);
    null->lc = null->rc = null;
}

void dfs(Node *u, bool init = true) {
    if (init) {
        u->f[0] = u;
        u->dep = 1;
        u->seg = null->insert(0, INT_MAX, u->val);
    }

    for (int i = 1; i < MAXN_LOG; i++) u->f[i] = u->f[i - 1]->f[i - 1];

    for (Edge *e = u->e; e; e = e->next) if (e->v != u->f[0]) {
        e->v->f[0] = u;
        e->v->dep = u->dep + 1;
        e->v->seg = u->seg->insert(0, INT_MAX, e->v->val);

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

int query(int u, int v, int k) {
    Node *p = lca(&N[u], &N[v]);
    PSegT *su = N[u].seg, *sv = N[v].seg;
    PSegT *sp = p->seg, *sf = (p == p->f[0] ? null : p->f[0]->seg);

    int l = 0, r = INT_MAX;
    while (l < r) {
        int mid = l + (r - l) / 2;
        int temp = su->lc->cnt + sv->lc->cnt - sp->lc->cnt - sf->lc->cnt;
        if (k <= temp) {
            su = su->lc;
            sv = sv->lc;
            sp = sp->lc;
            sf = sf->lc;
            r = mid;
        } else {
            k -= temp;
            su = su->rc;
            sv = sv->rc;
            sp = sp->rc;
            sf = sf->rc;
            l = mid + 1;
        }
    }
    return l;
}

int main() {
    init();

    int n, q;
    scanf("%d %d", &n, &q);

    for (int i = 1; i <= n; i++) scanf("%d", &N[i].val);

    for (int i = 1, u, v; i < n; i++) {
        scanf("%d %d", &u, &v);
        addEdge(u, v);
    }

    dfs(&N[1]);

    int lastAns = 0;
    while (q--) {
        int u, v, k;
        scanf("%d %d %d", &u, &v, &k);
        u ^= lastAns;
        printf("%d", lastAns = query(u, v, k));
        q ? puts("") : 0;
    }

    return 0;
}
```
