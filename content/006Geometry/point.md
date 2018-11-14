# 点相关

```c++
#include <cstdio>
#include <cmath>
#include <algorithm>

const int MAXN = 100005;
const double EPS = 1e-9;

int dcmp(double a, double b = 0) {
    double d = std::fabs(a - b);
    return d <= EPS ? 0 : (d > 0 ? 1 : -1);
}

struct Point {
    double x, y;
    
    Point(double x = 0, double y = 0) : x(x), y(y) {}
    
    bool operator<(const Point &rhs) const {
        return x == rhs.x ? y < rhs.y : x < rhs.x;
    }

    friend Point operator+(const Point &a, const Point &b) {
        return Point(a.x + b.x, a.y + b.y);
    }

    friend Point operator-(const Point &a, const Point &b) {
        return Point(a.x - b.x, a.y - b.y);
    }

    friend Point operator*(const Point &p, const double a) {
        return Point(p.x * a, p.y * a);
    }

    friend Point operator/(const Point &p, const double a) {
        return Point(p.x / a, p.y / a);
    }

    friend double dot(const Point &a, const Point &b) {
        return a.x * b.x + a.y * b.y;
    }

    friend double cross(const Point &a, const Point &b) {
        return a.x * b.y - a.y * b.x;
    }

    friend double angle(const Point &a, const Point &b) {
        return std::acos(dot(a, b) / a.length() / b.length());
    }

    Point rotate(double rad) {
        return Point(x * std::cos(rad) - y * std::sin(rad), x * std::sin(rad) + y * std::cos(rad));
    }

    double length() const {
        return std::sqrt(dot(*this, *this));
    }

    Point getPerpendicular() {
        double X = sqrt(1 / (1 + (x / y) * (x / y)));

        int sx, sy;
        if (x > 0 && y > 0) sx = 1, sy = -1;
        else if (x > 0 && y <= 0) sx = 1, sy = 1;
        else if (x <= 0 && y > 0) sx = -1, sy = -1;
        else sx = 1, sy = -1;

        return Point(sx * X, sy * sqrt(1 - X * X));
    }
} P[MAXN], ch[MAXN];

Point getLineIntersect(const Point &sa, const Point &ta, const Point &sb, const Point &tb) {
    double t = cross(tb - sb, sa - sb) / cross(ta - sa, tb - sb);
    return sa + (sb - sa) * t;
}

Point getLineProj(const Point &p, const Point &s, const Point &t) {
    Point u = t - s;
    return s + u * (dot(u, p - s) / dot(u, u));
}

double getDistToLine(const Point &p, const Point &s, const Point &t) {
    return std::fabs(cross(t - s, p - s)) / (t - s).length();
}

double getDistToSeg(const Point &p, const Point &s, const Point &t) {
    if (dcmp(dot(t - s, p - s)) < 0) return (p - s).length();
    else if (dcmp(dot(t - s, p - t)) > 0) return (p - t).length();
    else return std::fabs(cross(t - s, p - s)) / (t - s).length();
}

bool doesPointOnSeg(const Point &p, const Point &s, const Point &t) {
    return dcmp(cross(s - p, t - p)) == 0 && dcmp(dot(s - p, t - p)) < 0;
}

int getConvexHull(int n) {
    std::sort(P, P + n);

    int m = 0;
    for (int i = 0; i < n; i++) {
        while (m > 1 && cross(ch[m - 1] - ch[m - 2], P[i] - ch[m - 2]) <= 0) m--;
        ch[m++] = P[i];
    }

    int k = m;
    for (int i = n - 1; ~i; i--) {
        while (m > k && cross(ch[m - 1] - ch[m - 2], P[i] - ch[m - 2]) <= 0) m--;
        ch[m++] = P[i];
    }

    m > 1 ? m-- : 0;
    return m;
}

int getv(const Point &p){
    if (p.x >= 0 && p.y > 0) return 1;
    if (p.x < 0 && p.y >= 0) return 2;
    if (p.x <= 0 && p.y < 0) return 3;
    if (p.x > 0 && p.y <= 0) return 4;
}

bool cmpByAngle(const Point &a, const Point &b) {
    if (getv(a) == getv(b)) return cross(a, b) > 0;
    else return getv(a) < getv(b);
}

int main() {

    return 0;
}
```
