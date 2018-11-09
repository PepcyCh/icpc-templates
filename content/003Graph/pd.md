# 点分治

```c++
#include <cstdio>
#include <vector>
#include <stack>
#include <algorithm>

const int MAXN = 100005;

struct Node;
struct Edge;

struct Node {
    Edge *e;
    bool solved;
    int size, max;
} N[MAXN];

struct Edge {
    Node *u, *v;
    Edge *next;

    Edge() {}
    Edge(Node *u, Node *v) : u(u), v(v), next(u->e) {}
} _pool[MAXN << 1], *_curr = _pool;

void addEdge(int u, int v) {
    N[u].e = new (_curr++) Edge(&N[u], &N[v]);
    N[v].e = new (_curr++) Edge(&N[v], &N[u]);
}

void dfs(Node *u, Node *fa, std::vector<Node *> &vec) {
    u->size = 1;
    u->max = 0;
    vec.push_back(u);
    for (Edge *e = u->e; e; e = e->next) if (!e->v->solved && e->v != fa) {
        dfs(e->v, u, vec);
        u->size += e->v->size;
        u->max = std::max(u->max, e->v->size);
    }
}

Node *center(Node *s) {
    static std::vector<Node *> vec;
    dfs(s, NULL, vec);

    Node *res = NULL;
    for (int i = 0; i < vec.size(); i++) {
        vec[i]->max = std::max(vec[i]->max, s->size - vec[i]->size);
        if (!res || res->max > vec[i]->max) res = vec[i];
    }

    return res;
}

int calc(Node *u, Node *fa = NULL) {
    int res = 0;
    for (Edge *e = u->e; e; e = e->next) if (e->v != fa && !e->v->solved) {
        res += calc(e->v, u);
        // do something...
    }
    // do something...
    return res;
}

int solve() {
    static std::stack<Node *> s;
    s.push(&N[1]);

    int ans = 0;
    while (!s.empty()) {
        Node *u = s.top();
        s.pop();

        Node *root = center(u);
        root->solved = true;

        ans += calc(root);

        for (Edge *e = root->e; e; e = e->next) if (!e->v->solved) s.push(e->v);
    }

    return ans;
}

int main() {
    int n;
    scanf("%d", &n);

    for (int i = 1, u, v; i < n; i++) {
        scanf("%d %d", &u, &v);
        addEdge(u, v);
    }

    int ans = solve();
    printf("%d\n", ans);
}
```
