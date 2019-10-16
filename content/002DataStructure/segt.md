# 线段树（Segment Tree）

```c++
#include <cstdio>
#include <algorithm>

const int MAXN = 100005;

class SegT {
private:
    struct Node {
        Node *lc, *rc;
        int l, r;
        long long sum, tag;

        Node() {}
        Node(int pos, int val) : lc(NULL), rc(NULL), l(pos), r(pos), sum(val), tag(0) {}
        Node(Node *lc, Node *rc) : lc(lc), rc(rc), l(lc->l), r(rc->r), tag(0) {
            maintain();
        }

        void add(long long d) {
            sum += (r - l + 1) * d;
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
            sum = lc->sum + rc->sum;
        }

        void update(int l, int r, int d) {
            if (r < this->l || this->r < l) return;
            if (l <= this->l && this->r <= r) {
                add(d);
                return;
            }
            pushDown();
            lc->update(l, r, d);
            rc->update(l, r, d);
            maintain();
        }
        void update(int pos, int d) {
            if (l == r) {
                add(d);
                return;
            }
            pushDown();
            int mid = l + ((r - l) >> 1);
            (pos <= mid ? lc : rc)->update(pos, d);
            maintain();
        }

        long long query(int l, int r) {
            if (r < this->l || this->r < l) return 0;
            if (l <= this->l && this->r <= r) return sum;
            pushDown();
            return lc->query(l, r) + rc->query(l, r);
        }
        long long query(int pos) {
            if (l == r) return sum;
            pushDown();
            int mid = l + ((r - l) >> 1);
            return (pos <= mid ? lc : rc)->query(pos);
        }
    } *root, _pool[MAXN << 1], *_curr;

public:
    Node *build(int l, int r, int *a) {
        if (l == r) return new (_curr++) Node(l, a[l]);
        int mid = l + ((r - l) >> 1);
        return new (_curr++) Node(build(l, mid, a), build(mid + 1, r, a));
    }
    void build(int n, int *a) {
        _curr = _pool;
        root = build(1, n, a);
    }

    void update(int l, int r, int d) { root->update(l, r, d); }
    void update(int pos, int d) { root->update(pos, d); }
    long long query(int l, int r) { return root->query(l, r); }
    long long query(int pos) { return root->query(pos); }
} segT;

int main() {
    int n;
    scanf("%d", &n);

    static int a[MAXN];
    for (int i = 1; i <= n; i++) scanf("%d", &a[i]);

    segT.build(1, n, a);

    int q;
    scanf("%d", &q);
    while (q--) {
        int op, l, r;
        scanf("%d %d %d", &op, &l, &r);

        if (op == 1) {
            int d;
            scanf("%d", &d);
            segT.update(l, r, d);
        } else {
            printf("%lld\n", segT.query(l, r));
        }
    }
    
    return 0;
}
```
