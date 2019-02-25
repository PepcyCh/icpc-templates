# Deluanay 三角剖分与平面欧几里得距离最小生成树

```c++
#include <cstdio>
#include <cmath>
#include <set>
#include <list>
#include <vector>
#include <algorithm>

const int MAXN = 100005;
const double EPS = 1e-9;

int dcmp(double a, double b = 0) {
    double d = a - b;
    return std::abs(d) <= EPS ? 0 : (d > 0 ? 1 : -1);
}

double sqr(double x) { return x * x; }

struct Point {
    double x, y;
    int id;

    Point(double x = 0, double y = 0, int id = -1) : x(x), y(y), id(id) {}

    bool operator<(const Point &rhs) const {
        return dcmp(x, rhs.x) == 0 ? y < rhs.y : x < rhs.x;
    }

    Point operator+(const Point &rhs) const { return Point(x + rhs.x, y + rhs.y); }
    Point operator-(const Point &rhs) const { return Point(x - rhs.x, y - rhs.y); }
    friend double dot(const Point &a, const Point &b) { return a.x * b.x + a.y * b.y; }
    friend double cross(const Point &a, const Point &b) { return a.x * b.y - a.y * b.x; }

    friend double dist(const Point &a, const Point &b) {
        return std::sqrt(sqr(a.x - b.x) + sqr(a.y - b.y));
    }
    friend double distSqr(const Point &a, const Point &b) {
        return sqr(a.x - b.x) + sqr(a.y - b.y);
    }
} P[MAXN];

namespace Deluanay {

struct Point3D {
    double x, y, z;

    Point3D(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
    Point3D(const Point &p) : x(p.x), y(p.y), z(sqr(p.x) + sqr(p.y)) {}

    friend Point3D operator-(const Point3D &a, const Point3D &b) {
        return Point3D(a.x - b.x, a.y - b.y, a.z - b.z);
    }

    friend double dot(const Point3D &a, const Point3D &b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }
    friend Point3D cross(const Point3D &a, const Point3D &b) {
        return Point3D(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
    }
    friend double mix(const Point3D &a, const Point3D &b, const Point3D &c) {
        return dot(a, cross(b, c));
    }
};

struct Edge {
    int id;
    std::list<Edge>::iterator c;

    Edge(int id = 0) : id(id) {}
};

bool doesPointInCircle(const Point &a, Point b, Point c, const Point &p) {
    if (cross(b - a, c - a) < 0) std::swap(b, c);
    Point3D a3(a), b3(b), c3(c), p3(p);
    b3 = b3 - a3, c3 = c3 - a3, p3 = p3 - a3;
    return dcmp(mix(p3, b3, c3)) < 0;
}

bool hasIntersection(const Point &a, const Point &b, const Point &c, const Point &d) {
    return dcmp(cross(c - a, b - a)) * dcmp(cross(b - a, d - a)) > 0 &&
           dcmp(cross(a - c, d - c)) * dcmp(cross(d - c, b - c)) > 0;
}

std::list<Edge> head[MAXN];
Point p[MAXN];
int n, rename[MAXN];

void addEdge(int u, int v) {
    head[u].push_front(Edge(v));
    head[v].push_front(Edge(u));
    head[u].begin()->c = head[v].begin();
    head[v].begin()->c = head[u].begin();
}

void divide(int l, int r) {
    if (r - l <= 2) {
        for (int i = l; i <= r; i++) for (int j = i + 1; j <= r; j++) addEdge(i, j);
        return;
    }

    int mid = l + ((r - l) >> 1);
    divide(l, mid);
    divide(mid + 1, r);

    int nowl = l, nowr = r;
    for (bool updated = true; updated; ) {
        updated = false;
        Point ptL = p[nowl], ptR = p[nowr];
        for (auto i : head[nowl]) {
            Point t = p[i.id];
            double v = cross(ptL - ptR, t - ptR);
            if (dcmp(v) > 0 || (!dcmp(v) && distSqr(ptR, t) < distSqr(ptR, ptL))) {
                nowl = i.id;
                updated = true;
                break;
            }
        }
        if (updated) continue;
        for (auto i : head[nowr]) {
            Point t = p[i.id];
            double v = cross(ptR - ptL, t - ptL);
            if (dcmp(v) < 0 || (!dcmp(v) && distSqr(ptL, t) < distSqr(ptL, ptR))) {
                nowr = i.id;
                updated = true;
                break;
            }
        }
    }

    addEdge(nowl, nowr);
    while (true) {
        Point ptL = p[nowl], ptR = p[nowr];
        int ch = -1, side = 0;
        for (auto i : head[nowl]) {
            Point t = p[i.id];
            if (dcmp(cross(ptR - ptL, t - ptL)) > 0
                    && (ch == -1 || doesPointInCircle(ptL, ptR, p[ch], t))) {
                ch = i.id;
                side = -1;
            }
        }
        for (auto i : head[nowr]) {
            Point t = p[i.id];
            if (dcmp(cross(t - ptR, ptL - ptR)) > 0
                    && (ch == -1 || doesPointInCircle(ptL, ptR, p[ch], t))) {
                ch = i.id;
                side = 1;
            }
        }

        if (ch == -1) break;
        if (side == -1) {
            for (auto it = head[nowl].begin(); it != head[nowl].end(); ) {
                if (hasIntersection(ptL, p[it->id], ptR, p[ch])) {
                    head[it->id].erase(it->c);
                    head[nowl].erase(it++);
                } else {
                    it++;
                }
            }
            nowl = ch;
            addEdge(nowl, nowr);
        } else {
            for (auto it = head[nowr].begin(); it != head[nowr].end(); ) {
                if (hasIntersection(ptR, p[it->id], ptL, p[ch])) {
                    head[it->id].erase(it->c);
                    head[nowr].erase(it++);
                } else {
                    it++;
                }
            }
            nowr = ch;
            addEdge(nowl, nowr);
        }
    }
}

bool isParallel(const Point &a, const Point &b, const Point &c) {
    return !dcmp(cross(b - a, c - a));
}

void getEdge(std::vector<std::pair<int, int> > &ret) {
    ret.reserve(n);
    static std::set<std::pair<int, int> > vis;
    vis.clear();

    for (int i = 0; i < n; i++) for (auto j : head[i]) {
        if (j.id < i) continue;
        Point now = p[j.id];
        for (auto k : head[i]) {
            if (k.id < i) continue;
            if (isParallel(p[i], p[j.id], p[k.id])
                    && distSqr(p[k.id], p[i]) < distSqr(now, p[i]))
                now = p[k.id];
        }
        if (vis.find(std::make_pair(p[i].id, now.id)) == vis.end()) {
            vis.insert(std::make_pair(p[i].id, now.id));
            ret.push_back(std::make_pair(p[i].id, now.id));
        }
    }
}

void init(int n, Point *P) {
    std::copy(P, P + n, p);
    std::sort(p, p + n);
    for (int i = 0; i < n; i++) rename[p[i].id] = i;
    Deluanay::n = n;
    divide(0, n - 1);
}

void split(int n, Point *P, std::vector<std::pair<int, int> > &ret) {
    for (int i = 0; i < n; i++) P[i].id = i;
    init(n, P);
    getEdge(ret);
}

} // namespace Deluanay

namespace Kruskal {

struct DJS {
    int f[MAXN];

    void init(int n) { for (int i = 1; i <= n; i++) f[i] = i; }
    int find(int x) { return x == f[x] ? x : f[x] = find(f[x]); }
    void merge(int x, int y) { f[find(y)] = find(x); }
    bool test(int x, int y) { return find(x) == find(y); }
} djs;

struct Edge {
    int u, v;
    double w;

    Edge(int u, int v, double w) : u(u), v(v), w(w) {}
};
std::vector<Edge> edges;

double solve(int n) {
    std::sort(edges.begin(), edges.end(), [](const Edge &a, const Edge &b){
            return a.w < b.w;
        });
    djs.init(n);

    double ans = 0;
    int added = 0;
    for (auto e : edges) if (!djs.test(e.u, e.v)) {
        djs.merge(e.u, e.v);
        ans += e.w;
        if (++added == n - 1) break;
    }

    return ans;
}

} // namespace Kruskal

int main() {
    int n;
    scanf("%d", &n);
    for (int i = 0; i < n; i++) scanf("%lf %lf", &P[i].x, &P[i].y);

    static std::vector<std::pair<int, int> > edges;
    Deluanay::split(n, P, edges);

    for (auto e : edges) Kruskal::edges.emplace_back(e.first + 1, e.second + 1, dist(P[e.first], P[e.second]));
    double ans = Kruskal::solve(n);
    printf("%lld\n", ans);
    
    return 0;
}
```
