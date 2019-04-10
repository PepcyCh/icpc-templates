# 后缀自动机（Suffix Automaton）

```c++
#include <cstdio>
#include <vector>
#include <algorithm>

const int MAXN = 100005;
const int CHAR_SET = 26;

struct SAM {
    struct Node {
        Node *c[CHAR_SET], *next;
        int max, posCnt;

        Node(int max = 0, bool newSuffix = false) : max(max), posCnt(newSuffix), next(NULL), c() {}

        int min() const {
            return next->max + 1;
        }
    } *start, *last, _pool[MAXN << 1], *_curr;

    SAM() {
        init();
    }

    void init() {
        _curr = _pool;
        start = last = new (_curr++) Node;
    }

    Node *extend(int c) {
        Node *u = new (_curr++) Node(last->max + 1, true), *v = last;

        for (; v && !v->c[c]; v = v->next) v->c[c] = u;

        if (!v) {
            u->next = start;
        } else if (v->c[c]->max == v->max + 1) {
            u->next = v->c[c];
        } else {
            Node *n = new (_curr++) Node(v->max + 1), *o = v->c[c];
            std::copy(o->c, o->c + CHAR_SET, n->c);
            n->next = o->next;
            u->next = o->next = n;
            for (; v && v->c[c] == o; v = v->next) v->c[c] = n;
        }

        return last = u;
    }

    std::vector<Node *> topo;
    std::vector<Node *> toposort() {
        static int buc[MAXN << 1];
        int max = 0;
        for (Node *p = _pool; p != _curr; p++) {
            max = std::max(max, p->max);
            buc[p->max]++;
        }
        for (int i = 1; i <= max; i++) buc[i] += buc[i - 1];

        topo.resize(_curr - _pool);
        for (Node *p = _pool; p != _curr; p++) topo[--buc[p->max]] = p;

        return topo;
    }

    void calc() {
        toposort();

        for (int i = topo.size() - 1; i; i--) topo[i]->next->posCnt += topo[i]->posCnt;
    }
} sam;

int main() {
    
    return 0;
}
```
