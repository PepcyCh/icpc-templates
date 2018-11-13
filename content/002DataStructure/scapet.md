# 替罪羊树（Scapegoat Tree）

```c++
#include <cstdio>
#include <vector>
#include <algorithm>

const int MAXN = 100005;

template <typename T, size_t SIZE>
struct MemoryPool {
    char mem[sizeof (T) * SIZE], *top, *del[SIZE], **delTop;

    MemoryPool() {
        init();
    }

    void init() {
        top = mem;
        delTop = del;
    }

    void *alloc() {
        if (delTop != del) return *(--delTop);
        char *res = top;
        top += sizeof (T);
        return (void *) res;
    }

    void free(void *p) {
        *(delTop++) = (char *) p;
    }
};

struct ScapegoatTree {
    static constexpr double ALPHA = 0.7;

    struct Node {
        Node *c[2];
        int val, size, valid;
        bool isDeled;

        static MemoryPool<Node, MAXN> pool;

        Node() {}
        Node(int val) : val(val), size(1), valid(1), isDeled(false) {
            c[0] = c[1] = NULL;
        }

        void *operator new(size_t) {
            return pool.alloc();
        }
        void operator delete(void *p) {
            pool.free(p);
        }

        bool bad() {
            return (c[0] ? c[0]->size : 0) > ALPHA * size
                || (c[1] ? c[1]->size : 0) > ALPHA * size;
        }

        void maintain() {
            size = 1 + (c[0] ? c[0]->size : 0) + (c[1] ? c[1]->size : 0);
            valid = !isDeled + (c[0] ? c[0]->valid : 0) + (c[1] ? c[1]->valid : 0);
        }
    } *root;

    void dfs(Node *u, std::vector<Node *> &vec) {
        if (!u) return;
        dfs(u->c[0], vec);
        if (!u->isDeled) vec.push_back(u);
        dfs(u->c[1], vec);
        if (u->isDeled) delete u;
    }

    Node *build(std::vector<Node *> &vec, int l, int r) {
        if (l >= r) return NULL;
        int mid = l + ((r - l) >> 1);
        Node *o = vec[mid];
        o->c[0] = build(vec, l, mid);
        o->c[1] = build(vec, mid + 1, r);
        o->maintain();
        return o;
    }

    void rebuild(Node *&u) {
        static std::vector<Node *> vec;
        vec.clear();
        dfs(u, vec);
        u = build(vec, 0, vec.size());
    }

    void insert(int x, Node *&u) {
        if (!u) {
            u = new Node(x);
            return;
        }

        ++u->size;
        ++u->valid;
        if (x >= u->val) insert(x, u->c[1]);
        else insert(x, u->c[0]);
        if (u->bad()) rebuild(u);
    }
    void insert(int x) {
        insert(x, root);
    }

    void del(int rank, Node *u) {
        if (!u) return;
        if (!u->isDeled && rank == (u->c[0] ? u->c[0]->valid : 0) + 1) {
            u->isDeled = true;
            --u->valid;
            return;
        }
        --u->valid;
        if (rank <= (u->c[0] ? u->c[0]->valid : 0)) del(rank, u->c[0]);
        else del(rank - (u->c[0] ? u->c[0]->valid : 0) - !u->isDeled, u->c[1]);
    }
    void del(int rank) {
        del(rank, root);
    }

    int getRank(int x) {
        int res = 1;
        Node *u = root;
        while (u) {
            if (u->val >= x) u = u->c[0];
            else {
                res += (u->c[0] ? u->c[0]->valid : 0) + !u->isDeled;
                u = u->c[1];
            }
        }
        return res;
    }

    int findKth(int k) {
        Node *u = root;
        while (u) {
            if (!u->isDeled && (u->c[0] ? u->c[0]->valid : 0) + 1 == k) return u->val;
            if ((u->c[0] ? u->c[0]->valid : 0) >= k) u = u->c[0];
            else {
                k -= (u->c[0] ? u->c[0]->valid : 0) + !u->isDeled;
                u = u->c[1];
            }
        }
    }

    int pred(int x) {
        return findKth(getRank(x) - 1);
    }

    int succ(int x) {
        return findKth(getRank(x + 1));
    }
} scapeT;
MemoryPool<ScapegoatTree::Node, MAXN> ScapegoatTree::Node::pool;

int main() {
    int n;
    scanf("%d", &n);
    while (n--) {
        int op, x;
        scanf("%d %d", &op, &x);

        if (op == 1) scapeT.insert(x);
        if (op == 2) scapeT.del(scapeT.getRank(x));
        if (op == 3) printf("%d\n", scapeT.getRank(x));
        if (op == 4) printf("%d\n", scapeT.findKth(x));
        if (op == 5) printf("%d\n", scapeT.pred(x));
        if (op == 6) printf("%d\n", scapeT.succ(x));
    }

    return 0;
}
```
