# 欧拉回路（Eulerian path/circuit）

* DFS 版本（适用于有向图和无向图，代码为无向图）
* 无递归版本 - 有向图
* 无递归版本 - 无向图

## DFS 版本

```c++
#include <cstdio>
#include <vector>
#include <algorithm>

const int MAXN = 100005;
const int MAXM = 200005;

struct Edge;
struct Node {
    Edge *e;
    int deg;
} N[MAXN];

struct Edge {
    Node *u, *v;
    Edge *next, *rev;
    bool used;

    Edge() {}
    Edge(Node *u, Node *v) : u(u), v(v), next(u->e), used(false) {}
} _pool[MAXM << 1], *_curr = _pool;

void addEdge(int u, int v) {
    N[u].e = new (_curr++) Edge(&N[u], &N[v]);
    N[v].e = new (_curr++) Edge(&N[v], &N[u]);
    N[u].e->rev = N[v].e;
    N[v].e->rev = N[u].e;
    ++N[u].deg;
    ++N[v].deg;
}

namespace Euler {
    std::vector<Edge *> path;

    void dfs(Node *u) {
        for (Edge *e = u->e; e; e = e->next) if (!e->used) {
            e->used = e->rev->used = true;
            dfs(e->v);
            path.push_back(e);
        }
    }

    bool solve(int n) {
        int cnt = 0;
        for (int i = 1; i <= n; i++) cnt += N[i].deg % 2;
        if (cnt > 2) return false;
        int s = -1;
        if (cnt) for (int i = 1; i <= n; i++) if (N[i].deg % 2) {
            s = i;
            break;
        }
        s = (s == -1 ? 1 : s);
        dfs(&N[s]);
        std::reverse(path.begin(), path.end());
        return true;
    }
}

int main() {
    int n, m;
    scanf("%d %d", &n, &m);
    for (int i = 0, u, v; i < m; i++) {
        scanf("%d %d", &u, &v);
        addEdge(u, v);
    }

    if (Euler::solve(n)) {
        printf("%d", Euler::path[0]->u - N);
        for (int i = 0; i < Euler::path.size(); i++) {
            printf(" -> %d", Euler::path[i]->v - N);
        }
        puts("");
    } else {
        puts("-1");
    }

    return 0;
}
```

## 有向图欧拉回路

```c++
#include <cstdio>
#include <vector>
#include <stack>
#include <algorithm>

const int MAXN = 100005;
const int MAXM = 400005;

struct Edge;
struct Node {
    Edge *e, *curr;
    int ideg, odeg;
} N[MAXN];

struct Edge {
    Node *u, *v;
    Edge *next;

    Edge() {}
    Edge(Node *u, Node *v) : u(u), v(v), next(u->e) {}
} _pool[MAXM], *_curr = _pool;

void addEdge(int u, int v) {
    N[u].e = new (_curr++) Edge(&N[u], &N[v]);
    ++N[u].odeg;
    ++N[v].ideg;
}

namespace Hierholzer {
    std::stack<Node *> currPath;
    std::vector<Node *> euler;

    void hierholzer(int s) {
        currPath.push(&N[s]); // delete this line if a Eulerian Path is required
        Node *u = &N[s];
        do {
            if (u->odeg == 0) {
                euler.push_back(currPath.top());
                currPath.pop();
                if (!currPath.empty()) u = currPath.top();
            } else {
                Node *v = u->curr->v;
                u->curr = u->curr->next;
                --u->odeg;
                u = v;
                currPath.push(u);
            }
        } while (!currPath.empty());
    }

    bool solve(int n) {
        int cntS = 0, cntT = 0;
        for (int i = 1; i <= n; i++) {
            if (N[i].ideg + 1 == N[i].odeg) ++cntS;
            else if (N[i].odeg + 1 == N[i].ideg) ++cntT;
            else if (N[i].ideg != N[i].odeg) return false;
        }
        if (cntS > 1 || cntT > 1) return false;
        int s = -1;
        for (int i = 1; i <= n; i++) {
            if (N[i].ideg + 1 == N[i].odeg) s = i;
            N[i].curr = N[i].e;
        }
        s = (s == -1 ? 1 : s);
        hierholzer(s);
        return true;
    }
}

int main() {
    int n, m;
    scanf("%d %d", &n, &m);
    for (int i = 0, u, v; i < m; i++) {
        scanf("%d %d", &u, &v);
        addEdge(u, v);
    }

    if (Hierholzer::solve(n)) {
        for (int i = Hierholzer::euler.size() - 1; i >= 0; i--) {
            printf("%d", Hierholzer::euler[i] - N);
            if (i) printf(" -> ");
        }
        puts("");
    } else {
        puts("-1");
    }
    
    return 0;
}
```

## 无向图欧拉回路

```c++
#include <cstdio>
#include <list>
#include <vector>
#include <stack>
#include <algorithm>

const int MAXN = 100005;
const int MAXM = 400005;

struct Edge;
struct Node {
    std::list<Edge> e;
    int deg;
} N[MAXN];

struct Edge {
    Node *u, *v;
    std::list<Edge>::iterator rev;

    Edge() {}
    Edge(Node *u, Node *v) : u(u), v(v) {}
};

void addEdge(int u, int v) {
    N[u].e.emplace_front(&N[u], &N[v]);
    N[v].e.emplace_front(&N[v], &N[u]);
    N[u].e.front().rev = N[v].e.begin();
    N[v].e.front().rev = N[u].e.begin();
    ++N[u].deg;
    ++N[v].deg;
}

namespace Hierholzer {
    std::stack<Node *> currPath;
    std::vector<Node *> euler;

    void hierholzer(int s) {
        currPath.push(&N[s]); // delete this line if a Eulerian Path is required
        Node *u = &N[s];
        do {
            if (u->deg == 0) {
                euler.push_back(currPath.top());
                currPath.pop();
                if (!currPath.empty()) u = currPath.top();
            } else {
                Node *v = u->e.front().v;
                v->e.erase(u->e.front().rev);
                u->e.erase(u->e.begin());
                --u->deg;
                --v->deg;
                u = v;
                currPath.push(u);
            }
        } while (!currPath.empty());
    }

    bool solve(int n) {
        int cnt = 0;
        for (int i = 1; i <= n; i++) if (N[i].deg % 2)  ++cnt;
        if (cnt > 2) return false;
        int s = -1;
        if (cnt) for (int i = 1; i <= n; i++) if (N[i].deg % 2) {
            s = i;
            break;
        }
        s = (s == -1 ? 1 : s);
        hierholzer(s);
        return true;
    }
}

int main() {
    int n, m;
    scanf("%d %d", &n, &m);
    for (int i = 0, u, v; i < m; i++) {
        scanf("%d %d", &u, &v);
        addEdge(u, v);
    }

    if (Hierholzer::solve(n)) {
        for (int i = Hierholzer::euler.size() - 1; i >= 0; i--) {
            printf("%d", Hierholzer::euler[i] - N);
            if (i) printf(" -> ");
        }
        puts("");
    } else {
        puts("-1");
    }
    
    return 0;
}
```
