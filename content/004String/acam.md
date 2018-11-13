# AC 自动机（Aho-Corasick Automaton）

```c++
#include <cstdio>
#include <queue>
#include <algorithm>

const int MAXN = 100005;
const int CHAR_SET = 26;

struct ACAM {
    struct Node {
        Node *c[CHAR_SET], *next, *fail;
        bool isWord;
        int refCnt;

        Node(bool isWord = false) : c(), next(NULL), fail(NULL), isWord(isWord), refCnt(0) {}

        void apply() {
            refCnt++;
            if (next) next->apply();
        }
    } *root;

    ACAM() : root(new Node) {}

    Node *insert(char *begin, char *end) {
        Node **u = &root;
        for (char *p = begin; p != end; p++) {
            if (!(*u)) *u = new Node;
            u = &(*u)->c[*p - 'a'];
        }

        if (!(*u)) *u = new Node(true);
        else (*u)->isWord = true;

        return *u;
    }

    void build() {
        std::queue<Node *> q;
        q.push(root);
        root->fail = root;
        root->next = NULL;

        while (!q.empty()) {
            Node *u = q.front();
            q.pop();

            for (int i = 0; i < CHAR_SET; i++) {
                Node *&c = u->c[i];

                if (!c) {
                    c = u == root ? root : u->fail->c[i];
                    continue;
                }

                Node *v = u->fail;
                c->fail = u != root && v->c[i] ? v->c[i] : root;
                c->next = c->fail->isWord ? c->fail : c->fail->next;
                q.push(c);
            }
        }
    }

    void exec(char *l, char *r) {
        Node *u = root;
        for (char *p = l; p != r; p++) {
            u = u->c[*p - 'a'];
            if (u->isWord) u->apply();
            else if (u->next) u->next->apply();
        }
    }
} acam;

int main() {
    
    return 0;
}
```
