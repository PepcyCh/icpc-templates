# K-d Tree

求曼哈顿距离最近点对和欧几里得距离 $$k$$ 远点对。

### 曼哈顿距离最近点对

```c++
#include <cstdio>
#include <climits>
#include <algorithm>

const int MAXN = 500005;

struct Point {
    int x, y;
    Point(int x = 0, int y = 0) : x(x), y(y) {}
} P[MAXN];

int dist(const Point &a, const Point &b) {
    return abs(a.x - b.x) + abs(a.y - b.y);
}

struct KDTree {
    static bool cmp1(const Point &a, const Point &b) {
        return a.y < b.y ? a.y < b.y : (a.y == b.y && a.x < b.x);
    }
    static bool cmp2(const Point &a, const Point &b) {
        return a.x < b.x ? a.x < b.x : (a.x == b.x && a.y < b.y);
    }

    struct Node {
        Node *c[2];
        Point p, r1, r2;

        Node() {}
        Node (Point p) : p(p), r1(p), r2(p) {
            c[0] = c[1] = NULL;
        }

        void maintain() {
            if (c[0]) {
                r1.x = std::min(r1.x, c[0]->r1.x);
                r1.y = std::min(r1.y, c[0]->r1.y);
                r2.x = std::max(r2.x, c[0]->r2.x);
                r2.y = std::max(r2.y, c[0]->r2.y);
            }
            if (c[1]) {
                r1.x = std::min(r1.x, c[1]->r1.x);
                r1.y = std::min(r1.y, c[1]->r1.y);
                r2.x = std::max(r2.x, c[1]->r2.x);
                r2.y = std::max(r2.y, c[1]->r2.y);
            }
        }

        int dist(const Point &p) {
            int res = 0;
            if (p.x < r1.x || r2.x < p.x) res += p.x < r1.x ? r1.x - p.x : p.x - r2.x;
            if (p.y < r1.y || r2.y < p.y) res += p.y < r1.y ? r1.y - p.y : p.y - r2.y;
            return res;
        }

        void query(const Point &p, int &res) {
            res = std::min(res, ::dist(this->p, p));

            if (!(c[0] || c[1])) return;
            int k = c[0] && c[1] ? c[0]->dist(p) > c[1]->dist(p) : (c[0] ? 0 : 1);

            c[k]->query(p, res);
            if (c[k ^ 1] && c[k ^ 1]->dist(p) < res) c[k ^ 1]->query(p, res);
        }
    } *root, _pool[MAXN << 1], *_curr;

    KDTree() : root(NULL) {
        _curr = _pool;
    }

    Node *build(int l, int r, Point P[], int d = 0) {
        if (l > r) return NULL;
        if (l == r) return new (_curr++) Node(P[l]);

        int mid = l + ((r - l) >> 1);
        std::nth_element(P + l, P + mid, P + r + 1, d ? cmp1 : cmp2);

        Node *u = new (_curr++) Node(P[mid]);
        u->c[0] = build(l, mid - 1, P, d ^ 1);
        u->c[1] = build(mid + 1, r, P, d ^ 1);
        u->maintain();

        return u;
    }

    void insert(const Point &p) {
        Node **u = &root;
        int d = 0;
        while (*u) {
            int k = (d ? cmp1(p, (*u)->p) : cmp2(p, (*u)->p)) ^ 1;
            d ^= 1;
            (*u)->r1.x = std::min(p.x, (*u)->r1.x);
            (*u)->r1.y = std::min(p.y, (*u)->r1.y);
            (*u)->r2.x = std::max(p.x, (*u)->r2.x);
            (*u)->r2.y = std::max(p.y, (*u)->r2.y);
            u = &(*u)->c[k];
        }
        *u = new (_curr++) Node(p);
    }

    int query(const Point &p) {
        int res = INT_MAX;
        root->query(p, res);
        return res;
    }
} kdT;

int main() {
    int n, m;
    scanf("%d %d", &n, &m);
    for (int i = 1; i <= n; i++) scanf("%d %d", &P[i].x, &P[i].y);

    kdT.root = kdT.build(1, n, P);

    while (m--) {
        int op;
        Point p;
        scanf("%d %d %d", &op, &p.x, &p.y);
        if (op == 1) kdT.insert(p);
        else printf("%d\n", kdT.query(p));
    }

    return 0;
}
```

### 欧几里得距离 $$k$$ 远点对

```c++
#include <cstdio>
#include <climits>
#include <vector>
#include <queue>
#include <algorithm>

const int MAXN = 100005;

struct Point {
    int x, y;
    Point(int x = 0, int y = 0) : x(x), y(y) {}
} P[MAXN];

long long dist(const Point &a, const Point &b) {
    return (long long) (a.x - b.x) * (a.x - b.x) + (long long) (a.y - b.y) * (a.y - b.y);
}

struct KDTree {
    static bool cmp1(const Point &a, const Point &b) {
        return a.y < b.y || (a.y == b.y && a.x < b.x);
    }
    static bool cmp2(const Point &a, const Point &b) {
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    }

    struct Node {
        Node *c[2];
        Point p, r1, r2;

        Node() {}
        Node(Point p) : p(p), r1(p), r2(p) {
            c[0] = c[1] = NULL;
        }

        void maintain() {
            if (c[0]) {
                r1.x = std::min(r1.x, c[0]->r1.x);
                r1.y = std::min(r1.y, c[0]->r1.y);
                r2.x = std::max(r2.x, c[0]->r2.x);
                r2.y = std::max(r2.y, c[0]->r2.y);
            }
            if (c[1]) {
                r1.x = std::min(r1.x, c[1]->r1.x);
                r1.y = std::min(r1.y, c[1]->r1.y);
                r2.x = std::max(r2.x, c[1]->r2.x);
                r2.y = std::max(r2.y, c[1]->r2.y);
            }
        }

        long long dist(const Point &p) {
            return std::max(std::max(::dist(p, r1), ::dist(p, r2)),
                            std::max(::dist(p, Point(r1.x, r2.y)), ::dist(p, Point(r2.x, r1.y))));
        }

        void query(const Point &p, std::priority_queue<long long, std::vector<long long>, std::greater<long long> > &q) {
            long long d = ::dist(p, this->p);

            if (d > q.top()) q.pop(), q.push(d);
            if (!(c[0] || c[1])) return;

            long long dis[2] = {c[0] ? c[0]->dist(p) : INT_MIN,
                                c[1] ? c[1]->dist(p) : INT_MIN};
            int k = dis[0] < dis[1];

            c[k]->query(p, q);
            if (c[k ^ 1] && dis[k ^ 1] > q.top()) c[k ^ 1]->query(p, q);
        }
    } *root, _pool[MAXN], *_curr;

    KDTree() : root(NULL) {
        _curr = _pool;
    }

    Node *build(Point *l, Point *r, int d = 0) {
        if (l > r) return NULL;
        if (l == r) return new (_curr++) Node(*l);

        Point *mid = l + ((r - l) >> 1);
        std::nth_element(l, mid, r + 1, d ? cmp1 : cmp2);

        Node *u = new (_curr++) Node(*mid);
        u->c[0] = build(l, mid - 1, d ^ 1);
        u->c[1] = build(mid + 1, r, d ^ 1);
        u->maintain();

        return u;
    }

    long long query(Point P[], int n, int k) {
        std::priority_queue<long long, std::vector<long long>, std::greater<long long> > q;
        while (!q.empty()) q.pop();
        for (int i = 0; i < k << 1; i++) q.push(-1);
        for (int i = 1; i <= n; i++) root->query(P[i], q);
        return q.top();
    }
} kdT;

int main() {
    int n, k;
    scanf("%d %d", &n, &k);
    for (int i = 1; i <= n; i++) scanf("%d %d", &P[i].x, &P[i].y);

    kdT.root = kdT.build(P + 1, P + n);
    printf("%lld\n", kdT.query(P, n, k));

    return 0;
}
```
