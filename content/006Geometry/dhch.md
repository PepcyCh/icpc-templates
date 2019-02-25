# 动态半凸壳

Splay 维护动态上凸壳

```c++
#include <cstdio>
#include <algorithm>

const int MAXN = 1000005;
const double EPS = 1e-9;

int dcmp(double a, double b = 0) {
    double d = a - b;
    return std::abs(d) <= EPS ? 0 : (d > 0 ? 1 : -1);
}

struct Splay {
    struct Node {
        Node *c[2], *fa, *pred, *succ;
        double x, y;

        Node() {}
        Node(Node *fa, double x, double y) : x(x), y(y), c(), fa(fa), pred(NULL), succ(NULL) {}

        int relation() { return fa->c[1] == this; }
    } *root, _pool[MAXN], *_curr;

    Splay() : root(NULL), _curr(_pool) {}

    void rotate(Node *u) {
        Node *o = u->fa;
        int x = u->relation();

        u->fa = o->fa;
        if (u->fa) u->fa->c[o->relation()] = u;

        o->c[x] = u->c[x ^ 1];
        if (u->c[x ^ 1]) u->c[x ^ 1]->fa = o;

        u->c[x ^ 1] = o;
        o->fa = u;
    }

    Node *splay(Node *u, Node *targetFa = NULL) {
        while (u->fa != targetFa) {
            if (u->fa->fa == targetFa) rotate(u);
            else if (u->relation() == u->fa->relation()) rotate(u->fa), rotate(u);
            else rotate(u), rotate(u);
        }
        if (!targetFa) root = u;
        return u;
    }

    double predSlope(Node *u) {
        return u->pred ? (u->y - u->pred->y) / (u->x - u->pred->x) : 1.0 / 0.0;
    }

    double succSlope(Node *u) {
        return u->succ ? (u->y - u->succ->y) / (u->x - u->succ->x) : -1.0 / 0.0;
    }

    Node *insert(double x, double y) {
        if (!root) {
            root = new (_curr++) Node(NULL, x, y);
            return root;
        }

        Node **u = &root, *fa = NULL;

        while (*u && dcmp(x, (*u)->x)) {
            fa = *u;
            u = &(*u)->c[dcmp(x, (*u)->x) > 0];
        }

        if (*u) {
            if (dcmp((*u)->y, y) >= 0) return splay(*u);
            (*u)->y = y;
        } else {
            (*u) = new (_curr++) Node(fa, x, y);

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

        Node *v = *u;
        if (dcmp(predSlope(v), succSlope(v)) <= 0) {
            splay(v->pred);
            splay(v->succ, v->pred);
            v->pred->succ = v->succ;
            v->succ->pred = v->pred;
            v->succ->c[0] = NULL;
            return NULL;
        }

        while (v->pred && dcmp(predSlope(v->pred), predSlope(v)) <= 0)
            v->pred = v->pred->pred;
        if (v->pred) {
            splay(v->pred);
            splay(v, v->pred);
            v->pred->succ = v;
        }
        v->c[0] = NULL;

        while (v->succ && dcmp(succSlope(v->succ), succSlope(v)) >= 0)
            v->succ = v->succ->succ;
        if (v->succ) {
            splay(v->succ);
            splay(v, v->succ);
            v->succ->pred = v;
        }
        v->c[1] = NULL;

        return splay(v);
    }

    Node *find(double slope) {
        Node *u = root;
        while (true) {
            if (dcmp(predSlope(u), slope) < 0) u = u->c[0];
            else if (dcmp(succSlope(u), slope) > 0) u = u->c[1];
            else return u;
        }
    }
} splay;

int main() {
    
    return 0;
}
```
