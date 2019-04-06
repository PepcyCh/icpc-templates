# 平面图区域（Farmland）

```c++
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

const int MAXN = 200005;
const double EPS = 1e-9;
const double PI = std::acos(-1.0);

int dcmp(double a, double b = 0) {
    double d = a - b;
    return std::abs(d) <= EPS ? 0 : (d > 0 ? 1 : -1);
}

struct Point {
    int x, y;

    friend long long cross(const Point &a, const Point &b) {
        return 1ll * a.x * b.y - 1ll * a.y * b.x;
    }
} P[MAXN];

struct Edge {
    int u, v;
    double a;
    bool vis;

    bool operator<(const Edge &rhs) const { return a < rhs.a; }

    Edge(int u, int v, double a) : u(u), v(v), a(a), vis(false) {}
};
std::vector<Edge> N[MAXN];

void addEdge(int u, int v, double a1, double a2) {
    N[u].emplace_back(u, v, a1);
    N[v].emplace_back(v, u, a2);
}

int find(const std::vector<Edge> &e, double a) {
    double d = a + PI - EPS;
    if (dcmp(d, PI) > 0) d -= 2 * PI;
    int res = std::upper_bound(e.begin(), e.end(), d, [](double a, const Edge &b) {
                return a < b.a;
            }) - e.begin() - 1;
    if (res < 0) res += e.size();
    return res;
}

int main() {
    int n, m;
    scanf("%d %d", &n, &m);
    for (int i = 1; i <= n; i++) scanf("%d %d", &P[i].x, &P[i].y);
    for (int i = 0, u, v; i < m; i++) {
        scanf("%d %d", &u, &v);
        double a1 = std::atan2(P[v].y - P[u].y, P[v].x - P[u].x);
        double a2 = std::atan2(P[u].y - P[v].y, P[u].x - P[v].x);
        addEdge(u, v, a1, a2);
    }
    for (int i = 1; i <= n; i++) std::sort(N[i].begin(), N[i].end());

    static std::vector<long long> ans;
    for (int i = 1; i <= n; i++) {
        for (Edge &e : N[i]) if (!e.vis) {
            int u = i, v = e.v;
            long long S = cross(P[u], P[v]);
            double lasta = e.a;
            e.vis = true;
            do {
                int t = find(N[v], lasta);
                if (N[v][t].vis) {
                    S = 0;
                    break;
                }
                u = v;
                v = N[u][t].v;
                lasta = N[u][t].a;
                N[u][t].vis = true;
                S += cross(P[u], P[v]);
            } while (v != i);
            if (S > 0) ans.push_back(S);
        }
    }

    std::sort(ans.begin(), ans.end());
    printf("%zu\n", ans.size());
    for (auto i : ans) printf("%lld\n", i);

    return 0;
}
```
