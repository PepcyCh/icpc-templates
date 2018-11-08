# 可持久化 Treap（Persistent Treap）

```c++
#include <cstdio>
#include <cstdlib>
#include <algorithm>

const int MAXN = 100005;

struct PTreap {
    struct Node {
        Node *lc, *rc;
        int val, size;

        Node() {}
        Node(int val, Node *lc, Node *rc) : val(val), lc(lc), rc(rc), size((lc ? lc->size : 0) + 1 + (rc ? rc->size : 0)) {}
    } *root, _pool[MAXN * 20], *_curr;

    PTreap() : root(NULL), _curr(_pool) {}

    static int size(const Node *u) {
        return u ? u->size : 0;
    }

    Node *merge(Node *a, Node *b) {
        if (!a) return b;
        if (!b) return a;

        if (rand() % (size(a) + size(b)) < size(a)) {
            return new (_curr++) Node(a->val, a->lc, merge(a->rc, b));
        } else {
            return new (_curr++) Node(b->val, merge(a, b->lc), b->rc);
        }
    }

    std::pair<Node *, Node *> split(Node *u, int pos) {
        std::pair<Node *, Node *> res(NULL, NULL);
        if (!u) return res;

        if (pos <= size(u->lc)) {
            res = split(u->lc, pos);
            u = new (_curr++) Node(u->val, res.second, u->rc);
        } else {
            res = split(u->rc, pos - size(u->lc) - 1);
            u = new (_curr++) Node(u->val, u->lc, res.first);
        }

        return res;
    }
};

int main() {
    
    return 0;
}
```
