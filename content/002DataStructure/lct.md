# Link-Cut Tree

```c++
#include <cstdio>
#include <algorithm>

const int MAXN = 30005;

struct LinkCutTree {
    struct Node {
        Node *c[2], *fa, *pathFa;
        int val, max;
        bool rev;

        Node() {}
        Node(int val) : val(val), max(val), rev(false), fa(NULL), pathFa(NULL), c() {}

        int relation() { return fa->c[1] == this; }

        void maintain() {
            max = val;
            if (c[0]) max = std::max(max, c[0]->max);
            if (c[1]) max = std::max(max, c[1]->max);
        }

        void pushDown() {
            if (rev) {
                std::swap(c[0], c[1]);
                if (c[0]) c[0]->rev ^= 1;
                if (c[1]) c[1]->rev ^= 1;
                rev = false;
            }
        }

        void rotate() {
            std::swap(pathFa, fa->pathFa);

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
        }

        void splay() {
            while (fa) {
                if (fa->fa) fa->fa->pushDown();
                fa->pushDown();
                pushDown();

                if (!fa->fa) rotate();
                else if (fa->relation() == relation()) fa->rotate(), rotate();
                else rotate(), rotate();
            }
        }

        void evert() {
            access();
            splay();
            rev ^= 1;
        }

        void expose() {
            splay();
            pushDown();
            if (c[1]) {
                std::swap(c[1]->fa, c[1]->pathFa);
                c[1] = NULL;
                maintain();
            }
        }

        bool splice() {
            splay();
            if (!pathFa) return false;

            pathFa->expose();
            pathFa->c[1] = this;
            std::swap(fa, pathFa);
            fa->maintain();
            return true;
        }

        void access() {
            expose();
            while (splice());
        }

        int query() {
            access();
            splay();
            return max;
        }
    } *N[MAXN], _pool[MAXN], *_curr;

    LinkCutTree() : _curr(_pool) {}

    void makeTree(int u, int val) {
        N[u] = new (_curr++) Node(val);
    }

    void link(int u, int v) {
        N[v]->evert();
        N[v]->pathFa = N[u];
    }

    void cut(int u, int v) {
        N[u]->evert();
        N[v]->access();
        N[v]->splay();
        N[v]->pushDown();
        N[v]->c[0]->fa = NULL;
        N[v]->c[0] = NULL;
        N[v]->maintain();
    }

    int query(int u, int v) {
        N[u]->evert();
        return N[v]->query();
    }
} lct;

int main() {

    return 0;
}
```
