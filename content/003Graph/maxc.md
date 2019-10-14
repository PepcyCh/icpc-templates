# 最大团/最大独立集/弦图最小染色

最大团 = 补图最大独立集

对于弦图：最大团 = 弦图最小染色。

```c++
#include <cstdio>
#include <algorithm>

const int MAXN = 42;

bool G[MAXN][MAXN];

namespace MaxClique {

int n, ans, cnt[MAXN], group[MAXN], vis[MAXN];

bool dfs(int u, int pos) {
    for (int i = u + 1; i <= n; i++) {
        if (cnt[i] + pos <= ans) return false;
        if (G[u][i]) {
            int j;
            for (j = 0; j < pos; j++) if (!G[i][vis[j]]) break;
            if (j == pos) {
                vis[pos] = i;
                if (dfs(i, pos + 1)) return true;
            }
        }
    }

    if (pos > ans) {
        for (int i = 0; i < pos; i++) group[i] = vis[i];
        ans = pos;
        return true;
    }

    return false;
}

int solve(int n) {
    ans = -1;
    MaxClique::n = n;

    for (int i = n; i; i--) {
        vis[0] = i;
        dfs(i, 1);
        cnt[i] = ans;
    }

    return ans < 0 ? 0 : ans;
}

} // namespace MaxClique

int main() {
    int n, m;
    scanf("%d %d", &n, &m);

    for (int i = 1; i <= n; i++) std::fill(G[i] + 1, G[i] + n + 1, true);
    for (int i = 0, u, v; i < m; i++) {
        scanf("%d %d", &u, &v);
        G[u][v] = G[v][u] = false;
    }

    int ans = MaxClique::solve(n);
    printf("%d\n", ans);

    return 0;
}
```
