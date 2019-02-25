# 半平面交

```c++
#include <cstdio>
#include <cfloat>
#include <cmath>
#include <algorithm>

const int MAXN = 305;
const double EPS = 1e-14;

int dcmp(double a, double b = 0) {
    double d = a - b;
    return std::abs(d) <= EPS ? 0 : (d > 0 ? 1 : -1);
}

struct Point {
    double x, y;
    
    Point(double x = 0, double y = 0) : x(x), y(y) {}
    
    Point operator+(const Point &rhs) const { return Point(x + rhs.x, y + rhs.y); }
    Point operator-(const Point &rhs) const { return Point(x - rhs.x, y - rhs.y); }
    Point operator*(double rhs) const { return Point(x * rhs, y * rhs); }
    friend double cross(const Point &a, const Point &b) { return a.x * b.y - a.y * b.x; }
} P[MAXN], hpi[MAXN];

struct Line {
    Point p, v;
    double slop;
    
    Line() {}
    Line(const Point &p, const Point &v) : p(p), v(v) {
        slop = std::atan2(v.y, v.x);
    }
    
    Point getVal(double t) const { return p + v * t; }
    
    bool operator<(const Line &another) const {
        return slop < another.slop || (slop == another.slop && v.x
            && getVal(-p.x / v.x).y > getVal(-another.p.x / another.v.x).y);
    }
    
    friend Point getIntersection(const Line &a, const Line &b) {
        double t = cross(b.v, a.p - b.p) / cross(a.v, b.v);
        return a.getVal(t);
    }
} L[MAXN];

int n;
int halfplaneIntersection() {
    int cnt = 0;
    L[cnt++] = L[0];
    for (int i = 1; i <= n; i++) if (dcmp(L[i].slop, L[i - 1].slop)) L[cnt++] = L[i];
    std::sort(L, L + cnt);
    
    static Line q[MAXN];
    static Point p[MAXN];
    int l = 0, r = 0;
    q[l] = L[0];
    for (int i = 1; i < cnt; i++) {
        while (l < r && dcmp(cross(L[i].v, p[r - 1] - L[i].p)) < 0) r--;
        while (l < r && dcmp(cross(L[i].v, p[l] - L[i].p)) < 0) l++;
        q[++r] = L[i];
        if (l < r) p[r - 1] = getIntersection(q[r - 1], q[r]);
    }
    while (l < r && dcmp(cross(q[l].v, p[r - 1] - q[l].p)) < 0) r--;
    while (l < r && dcmp(cross(q[r].v, p[l] - q[r].p)) < 0) l++;
    
    if (r - l <= 1) return 0;
    
    cnt = 0;
    for (int i = l; i < r; i++) hpi[++cnt] = p[i];
    return cnt;
}

int main() {
    scanf("%d", &n);
    for (int i = 1; i <= n; i++) scanf("%lf", &P[i].x);
    for (int i = 1; i <= n; i++) scanf("%lf", &P[i].y);
    
    P[0] = Point(P[1].x, P[1].y + 1);
    P[n + 1] = Point(P[n].x, P[n].y + 1);
    for (int i = 0; i <= n; i++) L[i] = Line(P[i], P[i + 1] - P[i]);
    
    int m = halfplaneIntersection();
    
    return 0;
}
```
