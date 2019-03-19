# 二分图最大权完备匹配

Kuhn-Munkres 算法，复杂度 O(n^3)。

```c++
#include <cstdio>
#include <climits>
#include <algorithm>

const int MAXN = 505;

namespace KuhnMunkres {
    int n, G[MAXN][MAXN], x[MAXN], y[MAXN], prevx[MAXN], prevy[MAXN], mat[MAXN], slack[MAXN], par[MAXN];

    void adjust(int u) {
        mat[u] = prevy[u];
        if (prevy[mat[u]] != -2) adjust(prevx[mat[u]]);
    }

    bool find(int u) {
        for (int i = 1; i <= n; i++) {
            if (prevy[i] == -1) {
                if (slack[i] > x[u] + y[i] - G[u][i]) {
                    slack[i] = x[u] + y[i] - G[u][i];
                    par[i] = u;
                }
                if (x[u] + y[i] == G[u][i]) {
                    prevy[i] = u;
                    if (mat[i] == -1) {
                        adjust(i);
                        return true;
                    }
                    if (prevx[mat[i]] != -1) continue;
                    prevx[mat[i]] = i;
                    if (find(mat[i])) return true;
                }
            }
        }
        return false;
    }

    int solve(int n) {
        KuhnMunkres::n = n;
        std::fill(mat + 1, mat + n + 1, -1);
        std::fill(y + 1, y + n + 1, -1);
        for (int i = 1; i <= n; i++) {
            x[i] = 0;
            for (int j = 1; j <= n; j++) x[i] = std::max(x[i], G[i][j]);
        }
        bool flag = false;
        for (int i = 1; i <= n; i++) {
            std::fill(prevx + 1, prevx + n + 1, -1);
            std::fill(prevy + 1, prevy + n + 1, -1);
            std::fill(slack + 1, slack + n + 1, INT_MAX);
            prevx[i] = -2;
            if (find(i)) continue;
            flag = false;
            while (!flag) {
                int m = INT_MAX;
                for (int j = 1; j <= n; j++) if (prevy[j] == -1) m = std::min(m, slack[j]);
                for (int j = 1; j <= n; j++) {
                    if (prevx[j] != -1) x[j] -= m;
                    if (prevy[j] != -1) y[j] += m;
                    else slack[j] -= m;
                }
                for (int j = 1; j <= n; j++) {
                    if (prevy[j] == -1 && slack[j]) {
                        prevy[j] = par[j];
                        if (mat[j] == -1) {
                            adjust(j);
                            flag = true;
                            break;
                        }
                        prevx[mat[j]] = j;
                        if (find(mat[j])) {
                            flag = true;
                            break;
                        }
                    }
                }
            }
        }
        int res = 0;
        for (int i = 1; i <= n; i++) res += G[mat[i]][i];
        return res;
    }
} // namespace KuhnMunkres

int main() {
    int n1, n2, m;
    scanf("%d %d %d", &n1, &n2, &m);
    for (int i = 0, u, v, w; i < m; i++) {
        scanf("%d %d %d", &u, &v, &w);
        KuhnMunkres::G[u][v] = w;
    }
    printf("%d\n", KuhnMunkres::solve(std::max(n1, n2)));

    return 0;
}
```
