# 点相关

```c++
#include <cstdio>
#include <cmath>
#include <algorithm>

const int MAXN = 100005;
const double EPS = 1e-9;

int dcmp(double a, double b = 0) {
    double d = a - b;
    return std::abs(d) <= EPS ? 0 : (d > 0 ? 1 : -1);
}

struct Point {
    double x, y;
    
    Point(double x = 0, double y = 0) : x(x), y(y) {}
    
    bool operator<(const Point &rhs) const { return x == rhs.x ? y < rhs.y : x < rhs.x; }

    Point operator+(const Point &rhs) const { return Point(x + rhs.x, y + rhs.y); }
    Point operator-(const Point &rhs) const { return Point(x - rhs.x, y - rhs.y); }
    Point operator*(const double a) const { return Point(x * a, y * a); }
    Point operator/(const double a) const { return Point(x / a, y / a); }
    friend double dot(const Point &a, const Point &b) { return a.x * b.x + a.y * b.y; }
    friend double cross(const Point &a, const Point &b) { return a.x * b.y - a.y * b.x; }

    friend double angle(const Point &a, const Point &b) {
        return std::acos(dot(a, b) / a.length() / b.length());
    }

    Point rotate(double rad) {
        return Point(x * std::cos(rad) - y * std::sin(rad),
                     x * std::sin(rad) + y * std::cos(rad));
    }

    double length() const { return std::sqrt(dot(*this, *this)); }

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

bool parallel(const Point &as, const Point &at, const Point &bs, const Point &bt) {
    return dcmp(cross(at - as, bt - bs)) == 0;
}

bool getLineInter(const Point &as, const Point &at, const Point &bs, const Point &bt, Point &res) {
    if (parallel(as, at, bs, bt)) return false;
    double c1 = cross(as - bs, bt - bs);
    double c2 = cross(at - bs, bt - bs);
    res = (at * c1 - as * c2) / (c1 - c2);
    return true;
}

bool getSegInter(const Point &as, const Point &at, const Point &bs, const Point &bt, Point &p) {
	if (!dcmp(cross(at - as, bt - bs))) return false;
	double c1 = cross(bs - as, at - as), c2 = cross(bt - as, at - as);
	double c3 = cross(as - bs, bt - bs), c4 = cross(at - bs, bt - bs);
	if (dcmp(c1) * dcmp(c2) <= 0 && dcmp(c3) * dcmp(c4) <= 0) {
		p = (at * c3 - as * c4) / (c3 - c4);
		return true;
	} else return false;
}

Point getLineProj(const Point &p, const Point &s, const Point &t) {
    Point u = t - s;
    return s + u * (dot(u, p - s) / dot(u, u));
}

double getDistToLine(const Point &p, const Point &s, const Point &t) {
    return std::abs(cross(t - s, p - s)) / (t - s).length();
}

double getDistToSeg(const Point &p, const Point &s, const Point &t) {
    if (dcmp(dot(t - s, p - s)) < 0) return (p - s).length();
    else if (dcmp(dot(t - s, p - t)) > 0) return (p - t).length();
    else return std::abs(cross(t - s, p - s)) / (t - s).length();
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

// *r == *l should be satisfied for the following 3 functions
double polyArea(Point *l, Point *r) {
	double sum = 0;
	for (Point *p = l; p < r; p++)
        sum += cross(*(p + 1), *p);
	return std::abs(sum / 2.0);
}

bool insidePoly(Point *l, Point *r, const Point &p) {
	int num = 0;
	for (Point *q = l; q < r; q++) {
		if (doesPointOnSeg(p, *q, *(q + 1))) return true;
		int k = dcmp(cross(*(q + 1) - *q, p - *q));
		int d1 = dcmp(q->y - p.y);
		int d2 = dcmp((q + 1)->y - p.y);
		if (k > 0 && d1 <= 0 && d2 > 0) ++num;
		if (k < 0 && d2 <= 0 && d1 > 0) --num;
	}
	return num;
}

Point polyWeight(Point *l, Point *r, double area) {
	Point res(0, 0);
	for (Point *p = l; p < r; p++)
		res = res + (*p + *(p + 1)) * cross(*(p + 1), *p);
	return res / area / 6.0;
}

int main() {

    return 0;
}
```
