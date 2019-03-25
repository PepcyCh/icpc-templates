# 一般图最大匹配

带花树（Blossom Algorithm）。

```c++
#include <cstdio>
#include <vector>
#include <queue>
#include <algorithm>

const int MAXN = 505;

namespace Blossom {
    struct DJS{
        int f[MAXN];
        void init(int n) { for (int i = 1; i <= n; i++) f[i] = i; }
        int find(int x) { return x == f[x] ? x : f[x] = find(f[x]); }
        bool test(int x, int y) { return find(x) == find(y); }
    } djs;

	int n, match[MAXN], id[MAXN], fa[MAXN], vis[MAXN], tim;
	std::vector<int> G[MAXN];
    std::queue<int> q;

	int lca(int u, int v) {
        ++tim;
		while (vis[u] != tim) {
			if (u) {
				u = djs.find(u);
				if (vis[u] == tim) return u;
				vis[u] = tim;
				if (match[u]) u = djs.find(fa[match[u]]);
				else u = 0;
			}
			std::swap(u, v);
		}
		return u;
	}

	void update(int u, int v , int p) {
		while (djs.find(u) != p) {
			fa[u] = v;
			int t = match[u];
			if (id[t] == 1) id[t] = 0, q.push(t);
			if (djs.find(t) == t) djs.f[t] = p;
			if (djs.find(u) == u) djs.f[u] = p;
			v = t;
			u = fa[v];
		}
	}

	bool bfs(int s) {
		for (int i = 1; i <= n; i++) id[i] = -1;
        djs.init(n);
        while (!q.empty()) q.pop();

		q.push(s);
		id[s] = 0;
		while (!q.empty()) {
			int u = q.front();
			q.pop();
            for (int v : G[u]) {
				if (id[v] == -1) {
					fa[v] = u;
					id[v] = 1;
					if (!match[v]) {
						int curr = v;
						while (curr) {
							int temp = fa[curr], last = match[temp];
							match[temp] = curr;
							match[curr] = temp;
							curr = last;
						}
						return true;
					}
					id[match[v]] = 0;
					q.push(match[v]);
				} else if (id[v] == 0 && !djs.test(u, v)) {
					int p = lca(u, v);
					update(u, v, p);
					update(v, u, p);
				}
			}
		}
		return false;
	}

	int solve(int n) {
        Blossom::n = n;
        int res = 0;
		for (int i = 1; i <= n; i++) if (!match[i] && bfs(i)) ++res;
        return res;
	}
} // namespace Blossom

int main() {
    int n, m;
	scanf("%d %d", &n, &m);
	for (int i = 0, u, v; i < m; i++) {
		scanf("%d %d", &u, &v);
        Blossom::G[u].push_back(v);
        Blossom::G[v].push_back(u);
	}

	int ans = Blossom::solve(n);
    printf("%d\n", ans);
    for (int i = 1; i <= n; i++) printf("%d%c", Blossom::match[i], " \n"[i == n]);

	return 0;
}
```
