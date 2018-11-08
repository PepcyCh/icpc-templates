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
        Node *c[2], *fa, **root;
        int val, cnt, size;
        static MemoryPool<Node, MAXN> pool;

        Node() {}
        Node(Node **root, Node *fa, int val) : val(val), size(1), cnt(1), c(), fa(fa), root(root) {}

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

        void rotate() {
            Node *o = fa;
            int x = relation();

            fa = o->fa;
            if (fa) fa->c[o->relation()] = this;

            o->c[x] = c[x ^ 1];
            if (c[x ^ 1]) c[x ^ 1]->fa = o;

            c[x ^ 1] = o;
            o->fa = this;

            o->maintain();
            maintain();

            if (!fa) *root = this;
        }

        Node *splay(Node *targetFa = NULL) {
            while (fa != targetFa) {
                if (fa->fa == targetFa) rotate();
                else if (relation() == fa->relation()) fa->rotate(), rotate();
                else rotate(), rotate();
            }
            return this;
        }

        Node *pred() {
            Node *u = c[0];
            while (u->c[1]) u = u->c[1];
            return u;
        }

        Node *succ() {
            Node *u = c[1];
            while (u->c[0]) u = u->c[0];
            return u;
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
            (*u) = new Node(&root, fa, val);
        }

        return (*u)->splay();
    }

    void erase(Node *u) {
        Node *pred = u->pred(), *succ = u->succ();
        pred->splay();
        succ->splay(pred);

        if (u->cnt > 1) {
            --u->cnt;
            --u->size;
        } else {
            delete succ->c[0];
            succ->c[0] = NULL;
        }

        pred->size--;
        succ->size--;
    }
    void erase(int val) {
        Node *u = find(val);
        if (u) erase(u);
    }

    Node *find(int val) {
        Node *u = root;
        while (u && u->val != val) u = u->c[val > u->val];
        if (u) return u->splay();
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
        return u->splay();
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
            int res = u->pred()->val;
            erase(u);
            return res;
        }
        return u->pred()->val;
    }

    int succ(int val) {
        Node *u = find(val);
        if (!u) {
            u = insert(val);
            int res = u->succ()->val;
            erase(u);
            return res;
        }
        return u->succ()->val;
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
