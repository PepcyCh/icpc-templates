# 最小树形图

朱-刘算法。

```c++
#include <cstdio>
#include <climits>
#include <algorithm>

const int MAXN = 3005;

namespace ChuLiu {
    int G[MAXN][MAXN], n, eg[MAXN], queue[MAXN];
    bool used[MAXN], pass[MAXN], more;

    void combine(int id, long long &res) {
        int tot = 0, from;
        for (; id != 0 && !pass[id]; id = eg[id]) {
            queue[tot++] = id;
            pass[id] = true;
        }
        for (from = 0; from < tot && queue[from] != id; ++from) {}
        if (from == tot) return;
        more = true;
        for (int i = from; i < tot; i++) {
            res += G[eg[queue[i]]][queue[i]];
            if (i != from) {
                used[queue[i]] = true;
                for (int j = 1; j <= n; j++) if (!used[j])
                    G[id][j] = std::min(G[id][j], G[queue[i]][j]);
            }
        }
        for (int i = 1; i <= n; i++) if (!used[i] && i != id) {
            for (int j = from; j < tot; j++) {
                int k = queue[j];
                G[i][id] = std::min(G[i][id], G[i][k] - G[eg[k]][k]);
            }
        }
    }

    long long solve(int n, int root) {
        ChuLiu::n = n;
        std::fill(used + 1, used + n + 1, false);
        long long res = 0;
        for (more = true; more; ) {
            more = false;
            std::fill(eg + 1, eg + n + 1, 0);
            for (int i = 1, k; i <= n; i++) if (!used[i] && i != root) {
                k = 0;
                for (int j = 1; j <= n; j++) if (!used[j] && i != j) {
                    if (k == 0 || G[j][i] < G[k][i])
                        k = j;
                }
                eg[i] = k;
            }
            std::fill(pass + 1, pass + n + 1, false);
            for (int i = 1; i <= n; i++) if (!used[i] && !pass[i] && i != root)
                combine(i, res);
        }
        for (int i = 1; i <= n; i++) if (!used[i] && i != root) {
            if (G[eg[i]][i] == INT_MAX) return -1;
            res += G[eg[i]][i];
        }
        return res;
    }
}

int main() {
    int n, m, rt;
    scanf("%d %d %d", &n, &m, &rt);
    for (int i = 1; i <= n; i++)
        std::fill(ChuLiu::G[i] + 1, ChuLiu::G[i] + n + 1, INT_MAX);
    for (int i = 0, u, v, w; i < m; i++) {
        scanf("%d %d %d", &u, &v, &w);
        ChuLiu::G[u][v] = w;
    }

    long long ans = ChuLiu::solve(n, rt);
    printf("%lld\n", ans);

    return 0;
}
```
