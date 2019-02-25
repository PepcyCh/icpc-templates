# k 次圆覆盖

```c++
#include <cstdio>
#include <cmath>
#include <algorithm>

const int MAXN = 1005;
const double EPS = 1e-9;
const double PI = std::acos(-1.0);

int dcmp(double a, double b = 0.0) {
    double d = a - b;
    return std::abs(d) <= EPS ? 0 : (d > 0 ? 1 : -1);
}

inline double sqr(double x) { return x * x; }
inline double safeSqrt(double x) { return !dcmp(x) ? 0 : std::sqrt(x); }

struct Point {
    double x, y;

    Point(double x = 0, double y = 0) : x(x), y(y) {}

    Point operator+(const Point &rhs) const { return Point(x + rhs.x, y + rhs.y); }
    Point operator-(const Point &rhs) const { return Point(x - rhs.x, y - rhs.y); }
    friend double cross(const Point &a, const Point &b) { return a.x * b.y - a.y * b.x; }
    double length() const { return std::sqrt(x * x + y * y); }
};

double dist(const Point &a, const Point &b) {
    return std::sqrt(sqr(a.x - b.x) + sqr(a.y - b.y));
}

struct Circle {
    Point p;
    double r, angle;
    int d;

    Circle() {}
    Circle(const Point &p, double r) : p(p), r(r), d(1), angle(0) {}
    Circle(double x, double y, double angle = 0, int d = 0) : p(Point(x, y)), r(0), d(d), angle(angle) {}
} C[MAXN];

namespace CircleUnion {

Circle tp[MAXN << 1];
double area[MAXN];

double calc(const Circle &c, const Circle &cp1, const Circle &cp2) {
    double ans = (cp2.angle - cp1.angle) * sqr(c.r) - cross(cp1.p - c.p, cp2.p - c.p) + cross(cp1.p, cp2.p);
    return ans / 2.0;
}

int numberOfCircleCross(const Circle a, const Circle b, Circle &cp1, Circle &cp2) {
    double mx = b.p.x - a.p.x, sx = b.p.x + a.p.x, mx2 = mx * mx;
    double my = b.p.y - a.p.y, sy = b.p.y + a.p.y, my2 = my * my;
    double sq = mx2 + my2, d = -(sq - sqr(a.r - b.r)) * (sq - sqr(a.r + b.r));

    if (dcmp(d) < 0) return 0;
    d = safeSqrt(d);

    double x = mx * ((a.r + b.r) * (a.r - b.r) + mx * sx) + sx * my2;
    double y = my * ((a.r + b.r) * (a.r - b.r) + my * sy) + sy * mx2;
    double dx = mx * d, dy = my * d;
    sq *= 2;

    cp1.p.x = (x - dy) / sq;
    cp1.p.y = (y + dx) / sq;
    cp2.p.x = (x + dy) / sq;
    cp2.p.y = (y - dx) / sq;
    return dcmp(d) > 0 ? 2 : 1;
}

void circleUnion(Circle *C, int n) {
    std::sort(C, C + n, [](const Circle &a, const Circle &b) {return a.r < b.r;});

    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++)
        if (dcmp(dist(C[i].p, C[j].p) + C[i].r - C[j].r) <= 0) ++C[i].d;

    for (int i = 0; i < n; i++) {
        int tn = 0, cnt = 0;
        for (int j = 0; j < n; j++) if (i != j) {
            static Circle cp1, cp2;
            // pay attention to the order of parameter!!!
            if (numberOfCircleCross(C[i], C[j], cp2, cp1) < 2) continue;
            cp1.angle = std::atan2(cp1.p.y - C[i].p.y, cp1.p.x - C[i].p.x);
            cp2.angle = std::atan2(cp2.p.y - C[i].p.y, cp2.p.x - C[i].p.x);
            cp1.d = 1;
            cp2.d = -1;
            tp[tn++] = cp1;;
            tp[tn++] = cp2;
            if (dcmp(cp1.angle, cp2.angle) > 0) ++cnt;
        }
        tp[tn++] = Circle(C[i].p.x - C[i].r, C[i].p.y, PI, -cnt);
        tp[tn++] = Circle(C[i].p.x - C[i].r, C[i].p.y, -PI, cnt);
        std::sort(tp, tp + tn, [](const Circle &a, const Circle &b) {
            return !dcmp(a.angle, b.angle) ? a.d > b.d : a.angle < b.angle;
        });
        int s = C[i].d + tp[0].d;
        for (int j = 1; j < tn; j++) {
            int p = s;
            s += tp[j].d;
            area[p] += calc(C[i], tp[j - 1], tp[j]);
        }
    }
}

void solve(Circle *C, int n) {
    std::fill(area + 1, area + n + 1, 0);
    circleUnion(C, n);

    // delete this line if you want to know area coverd by at least i-times instead of exactly i-times
    for (int i = 1; i < n; i++) area[i] -= area[i + 1];
}

} // namespace CircleUnion

int main() {
    int n;
    scanf("%d", &n);
    for (int i = 0; i < n; i++) {
        scanf("%lf %lf %lf", &C[i].p.x, &C[i].p.y, &C[i].r);
        C[i].d = 1;
    }

    CircleUnion::solve(C, n);
    
    return 0;
}
```
