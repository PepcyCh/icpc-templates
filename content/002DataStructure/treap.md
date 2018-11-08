# 树堆（Treap）

### 无旋式 Treap

```c++
#include <cstdio>
#include <cstdlib>
#include <algorithm>

const int MAXN = 100005;

struct Treap {
    struct Node {
        Node *lc, *rc;
        int val, size;
        bool rev;

        Node() {}
        Node(int val, Node *lc = NULL, Node *rc = NULL) : val(val), lc(lc), rc(rc), rev(false) {
            maintain();
        }

        void reverse() {
            std::swap(lc, rc);
            rev ^= 1;
        }

        void pushDown() {
            if (rev) {
                if (lc) lc->reverse();
                if (rc) rc->reverse();
                rev = false;
            }
        }

        void maintain() {
            size = (lc ? lc->size : 0) + 1 + (rc ? rc->size : 0);
        }

        void print() {
            pushDown();
            if (lc) lc->print();
            printf("%d ", val);
            if (rc) rc->print();
        }
    } *root, _pool[MAXN], *_curr;

    Treap() : root(NULL), _curr(_pool) {}

    static int size(const Node *u) {
        return u ? u->size : 0;
    }

    Node *_build(int l, int r) {
        if (l > r) return NULL;
        int mid = l + ((r - l) >> 1);
        return new (_curr++) Node(mid, _build(l, mid - 1), _build(mid + 1, r));
    }
    void build(int l, int r) {
        root = _build(l, r);
    }

    Node *merge(Node *a, Node *b) {
        if (!a) return b;
        if (!b) return a;

        if (rand() % (a->size + b->size) < a->size) {
            a->pushDown();
            a->rc = merge(a->rc, b);
            a->maintain();
            return a;
        } else {
            b->pushDown();
            b->lc = merge(a, b->lc);
            b->maintain();
            return b;
        }
    }

    std::pair<Node *, Node *> split(Node *u, int pos) {
        std::pair<Node *, Node *> res(NULL, NULL);
        if (!u) return res;
        u->pushDown();
        if (pos <= size(u->lc)) {
            res = split(u->lc, pos);
            u->lc = res.second;
            u->maintain();
            res.second = u;
        } else {
            res = split(u->rc, pos - size(u->lc) - 1);
            u->rc = res.first;
            u->maintain();
            res.first = u;
        }
        return res;
    }

    void reverse(int l, int r) {
        std::pair<Node *, Node *> L = split(root, l - 1);
        std::pair<Node *, Node *> R = split(L.second, r - l + 1);
        R.first->reverse();
        root = merge(merge(L.first, R.first), R.second);
    }

    void print() {
        root->print();
        puts("");
    }
} treap;

int main() {
    int n, m;
    scanf("%d %d", &n, &m);

    treap.build(1, n);

    for (int i = 0, l, r; i < m; i++) {
        scanf("%d %d", &l, &r);
        treap.reverse(l, r);
    }

    treap.print();

    return 0;
}
```
