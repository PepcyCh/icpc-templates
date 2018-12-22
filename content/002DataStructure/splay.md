# 伸展树（Splay）

```c++
#include <cstdio>
#include <climits>
#include <algorithm>

const int MAXN = 100005;

template <typename T, size_t SIZE>
struct MemoryPool {
    char mem[sizeof (T) * SIZE], *memTop, *del[SIZE], **delTop;

    MemoryPool() : memTop(mem), delTop(del) {}

    void *alloc() {
        if (delTop != del) return (void *) *--delTop;
        char *res = memTop;
        memTop += sizeof (T);
        return (void *) res;
    }

    void free(void *p) {
        *delTop++ = (char *) p;
    }
};

struct Splay {
    struct Node {
        Node *c[2], *fa, *pred, *succ;
        int val, cnt, size;
        static MemoryPool<Node, MAXN> pool;

        Node() {}
        Node(Node *fa, int val) : val(val), size(1), cnt(1), c(), fa(fa), pred(NULL), succ(NULL) {}

        ~Node() {
            if (c[0]) delete c[0];
            if (c[1]) delete c[1];
        }

        void *operator new(size_t) {
            return pool.alloc();
        }

        void operator delete(void *p) {
            pool.free(p);
        }

        void maintain() {
            size = (c[0] ? c[0]->size : 0) + cnt + (c[1] ? c[1]->size : 0);
        }

        int relation() {
            return fa->c[1] == this;
        }
    } *root;

    Splay() : root(NULL) {}

    static int size(const Node *u) {
        return u ? u->size : 0;
    }

    void init() {
        insert(INT_MIN);
        insert(INT_MAX);
    }

    void rotate(Node *u) {
        Node *o = u->fa;
        int x = u->relation();

        u->fa = o->fa;
        if (u->fa) u->fa->c[o->relation()] = u;

        o->c[x] = u->c[x ^ 1];
        if (u->c[x ^ 1]) u->c[x ^ 1]->fa = o;

        u->c[x ^ 1] = o;
        o->fa = u;

        o->maintain();
        u->maintain();

        if (!u->fa) root = u;
    }

    Node *splay(Node *u, Node *targetFa = NULL) {
        while (u->fa != targetFa) {
            if (u->fa->fa == targetFa) rotate(u);
            else if (u->relation() == u->fa->relation()) rotate(u->fa), rotate(u);
            else rotate(u), rotate(u);
        }
        return u;
    }

    Node *insert(int val) {
        Node **u = &root, *fa = NULL;

        while (*u && (*u)->val != val) {
            fa = *u;
            fa->size++;
            u = &(*u)->c[val > (*u)->val];
        }

        if (*u) {
            ++(*u)->cnt;
            ++(*u)->size;
        } else {
            (*u) = new Node(fa, val);

            if (fa) {
                if ((*u)->relation()) {
                    (*u)->succ = fa->succ;
                    (*u)->pred = fa;
                    if (fa->succ) fa->succ->pred = *u;
                    fa->succ = *u;
                } else {
                    (*u)->pred = fa->pred;
                    (*u)->succ = fa;
                    if (fa->pred) fa->pred->succ = *u;
                    fa->pred = *u;
                }
            }
        }

        return splay(*u);
    }

    void erase(Node *u) {
        splay(u->pred);
        splay(u->succ, u->pred);

        if (u->cnt > 1) {
            --u->cnt;
            --u->size;
        } else {
            delete u->succ->c[0];
            u->succ->c[0] = NULL;
            u->pred->succ = u->succ;
            u->succ->pred = u->pred;
        }

        --u->pred->size;
        --u->succ->size;
    }
    void erase(int val) {
        Node *u = find(val);
        if (u) erase(u);
    }

    Node *find(int val) {
        Node *u = root;
        while (u && u->val != val) u = u->c[val > u->val];
        if (u) return splay(u);
        else return NULL;
    }

    Node *select(int k) {
        Node *u = root;
        while (k < size(u->c[0]) || k >= size(u->c[0]) + u->cnt) {
            if (k < size(u->c[0])) u = u->c[0];
            else {
                k -= size(u->c[0]) + u->cnt;
                u = u->c[1];
            }
        }
        return splay(u);
    }

    int rank(int val) {
        Node *u = find(val);
        if (!u) {
            u = insert(val);
            int res = size(u->c[0]);
            erase(u);
            return res;
        }
        return size(u->c[0]);
    }

    int pred(int val) {
        Node *u = find(val);
        if (!u) {
            u = insert(val);
            int res = u->pred->val;
            erase(u);
            return res;
        }
        return u->pred->val;
    }

    int succ(int val) {
        Node *u = find(val);
        if (!u) {
            u = insert(val);
            int res = u->succ->val;
            erase(u);
            return res;
        }
        return u->succ->val;
    }
} splay;

MemoryPool<Splay::Node, MAXN> Splay::Node::pool;

int main() {
    splay.init();

    int n;
    scanf("%d", &n);

    while (n--) {
        int op, x;
        scanf("%d %d", &op, &x);

        if (op == 1) splay.insert(x);
        if (op == 2) splay.erase(x);
        if (op == 3) printf("%d\n", splay.rank(x));
        if (op == 4) printf("%d\n", splay.select(x)->val);
        if (op == 5) printf("%d\n", splay.pred(x));
        if (op == 6) printf("%d\n", splay.succ(x));
    }

    return 0;
}
```
