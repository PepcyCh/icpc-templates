# 二分图最大匹配

Hopcroft-Karp 算法，复杂度 O(E \sqrt{V})。

```c++
#include <cstdio>
#include <vector>
#include <queue>
#include <algorithm>

const int MAXN = 1005; // one side
const int MAXM = 1005; // the other side

namespace HopcroftKarp {
    std::vector<int> G[MAXN];
    int matx[MAXN], maty[MAXM], dx[MAXN], dy[MAXM];
    bool vis[MAXM];

    bool find(int u) {
        for (int v : G[u]) if (!vis[v] && dy[v] == dx[u] + 1) {
            vis[v] = true;
            if (!maty[v] || find(maty[v])) {
                matx[u] = v;
                maty[v] = u;
                return true;
            }
        }
        return false;
    }

    int solve(int n, int m) {
        std::fill(matx + 1, matx + n + 1, 0);
        std::fill(maty + 1, maty + m + 1, 0);
        int res = 0;
        while (true) {
            static std::queue<int> q;
            while (!q.empty()) q.pop();
            bool flag = false;
            std::fill(dx + 1, dx + n + 1, 0);
            std::fill(dy + 1, dy + m + 1, 0);
            for (int i = 1; i <= n; i++) if (!matx[i]) q.push(i);
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                for (int v : G[u]) if (!dy[v]) {
                    dy[v] = dx[u] + 1;
                    if (maty[v]) {
                        dx[maty[v]] = dy[v] + 1;
                        q.push(maty[v]);
                    } else {
                        flag = true;
                    }
                }
            }
            if (!flag) break;
            std::fill(vis + 1, vis + m + 1, false);
            for (int i = 1; i <= n; i++) if (!matx[i] && find(i)) ++res;
        }
        return res;
    }
} // namespace HopcroftKarp

int main() {
    int n1, n2, m;
    scanf("%d %d %d", &n1, &n2, &m);
    for (int i = 0, u, v; i < m; i++) {
        scanf("%d %d", &u, &v);
        HopcroftKarp::G[u].push_back(v);
    }
    printf("%d\n", HopcroftKarp::solve(n1, n2));
    
    return 0;
}
```
