# 三维凸包

```c++
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>

const int MAXN = 505;
const double EPS = 1e-9;

int dcmp(double a, double b = 0) {
    double d = a - b;
    return std::abs(d) <= EPS ? 0 : (d > 0 ? 1 : -1);
}

double asqrt(double x) { return x <= EPS ? 0 : std::sqrt(x); }
double sqr(double x) { return x * x; }

struct Point {
    double x, y, z;

    Point(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

    Point operator+(const Point &rhs) const {
        return Point(x + rhs.x, y + rhs.y, z + rhs.z);
    }
    Point operator-(const Point &rhs) const {
        return Point(x - rhs.x, y - rhs.y, z - rhs.z);
    }
    bool operator==(const Point &rhs) const {
        return !dcmp(x, rhs.x) && !dcmp(y, rhs.y) && !dcmp(z, rhs.z);
    }
    bool operator<(const Point &rhs) const {
        if (dcmp(x, rhs.x)) return x < rhs.x;
        if (dcmp(y, rhs.y)) return y < rhs.y;
        return z < rhs.z;
    }

    friend double dot(const Point &a, const Point &b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }
    friend Point cross(const Point &a, const Point &b) {
        return Point(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
    }
    friend double mix(const Point &a, const Point &b, const Point &c) {
        return dot(a, cross(b, c));
    }

    double length() const {
        return std::sqrt(sqr(x) + sqr(y) + sqr(z));
    }
} P[MAXN];

double volumn(const Point &a, const Point &b, const Point &c, const Point &d) {
    return mix(b - a, c - a, d - a) / 6.0;
}

bool doesPointOnLine(const Point &p, const Point &s, const Point &t) {
    return dcmp(cross(p - s, t - s).length()) == 0;
}

namespace ConvexHull3D {

struct Face {
    int a, b, c;
    bool judged;

    Face(int a = 0, int b = 0, int c = 0) : a(a), b(b), c(c), judged(false) {}

    int &operator[](int i) {
        return i == 0 ? a : (i == 1 ? b : c);
    }
    const int &operator[](int i) const {
        return i == 0 ? a : (i == 1 ? b : c);
    }
};
std::vector<Face> ch;

bool same(int a, int b) {
    return !dcmp(volumn(P[ch[a][0]], P[ch[a][1]], P[ch[a][2]], P[ch[b][0]]))
        && !dcmp(volumn(P[ch[a][0]], P[ch[a][1]], P[ch[a][2]], P[ch[b][1]]))
        && !dcmp(volumn(P[ch[a][0]], P[ch[a][1]], P[ch[a][2]], P[ch[b][2]]));
}

bool doesFaceHaveEdge(const Face &f, int s, int t) {
    for (int i = 0; i < 3; i++) {
        if (f[i] == s && f[(i + 1) % 3] == t) return true;
        if (f[i] == t && f[(i + 1) % 3] == s) return true;
    }
    return false;
}

bool find(int n) {
    for (int i = 2; i < n; i++) {
        if (doesPointOnLine(P[i], P[0], P[1])) continue;
        std::swap(P[2], P[i]);
        for (int j = i + 1; j < n; j++) {
            if (dcmp(volumn(P[0], P[1], P[2], P[j]))) {
                std::swap(P[3], P[j]);
                ch.emplace_back(0, 1, 2);
                ch.emplace_back(0, 2, 1);
                return true;
            }
        }
    }
    return false;
}

int mark[MAXN][MAXN];

void add(int pi, int &cnt) {
    static std::vector<Face> temp;
    temp.clear();

    ++cnt;
    for (auto f : ch) {
        int a = f[0], b = f[1], c = f[2];
        if (dcmp(volumn(P[pi], P[a], P[b], P[c])) < 0) {
            mark[a][b] = mark[b][a] = cnt;
            mark[a][c] = mark[c][a] = cnt;
            mark[b][c] = mark[c][b] = cnt;
        } else {
            temp.push_back(f);
        }
    }
    ch = temp;
    for (auto f : temp) {
        int a = f[0], b = f[1], c = f[2];
        if (mark[a][b] == cnt) ch.emplace_back(b, a, pi);
        if (mark[b][c] == cnt) ch.emplace_back(c, b, pi);
        if (mark[c][a] == cnt) ch.emplace_back(a, c, pi);
    }
}

bool getConvexHull(int n) {
    std::sort(P, P + n);
    n = std::unique(P, P + n) - P;
    std::random_shuffle(P, P + n);

    if (find(n)) {
        memset(mark, 0, sizeof (mark));
        int cnt = 0;
        for (int i = 3; i < n; i++) add(i, cnt);
    } else {
        return false;
    }

    // E, F, V: The number of edges, faces and vertices when each face is a triangle
    // EE, FF, VV: The number of edges, faces and vertices after dealing with co-planar triangles
    int F = ch.size(), V = (4 + F) / 2, E = V + F - 2;
    int FF = 0, EE = E;
    for (int i = 0; i < ch.size(); i++) if (!ch[i].judged) {
        static std::vector<Face> f;
        f.clear();
        for (int j = 0; j < ch.size(); j++) if (same(i, j) && !ch[j].judged) {
            f.push_back(ch[j]);
            ch[j].judged = true;
        }
        for (int j = 0; j < f.size(); j++) for (int k = 0; k < 3; k++) for (int l = 0; l < j; l++)
            if (doesFaceHaveEdge(f[l], f[j][k], f[j][(k + 1) % 3])) {
                --EE;
                break;
            }
        ++FF;
    }
    int VV = EE + 2 - FF;

    double area = 0; // surface area
    for (auto f : ch) {
        Point p = cross(P[f[0]] - P[f[1]], P[f[2]] - P[f[1]]);
        area += p.length() / 2.0;
    }

    double vol = 0; // volumn
    const Point &vp = P[ch[0][0]];
    for (auto f : ch) vol += std::abs(volumn(vp, P[f[0]], P[f[1]], P[f[2]]));

    return true;
}

} // namespace ConvexHull3D

int main() {
    int n;
    scanf("%d", &n);
    for (int i = 0; i < n; i++) scanf("%lf %lf %lf", &P[i].x, &P[i].y, &P[i].z);

    // Deal with co-linear points
    static bool banned[MAXN];
    for (int i = 0; i < n; i++) if (!banned[i]) for (int j = 0; j < i; j++) if (!banned[j]) {
        static std::vector<std::pair<Point, int> > l;
        l.clear();
        for (int k = 0; k < n; k++) if (!banned[k] && doesPointOnLine(P[k], P[i], P[j]))
            l.emplace_back(P[k], k);
        std::sort(l.begin(), l.end());
        for (int k = 1; k < l.size() - 1; k++) banned[l[k].second] = true;
    }

    int p = 0;
    for (int i = 0; i < n; i++) if (!banned[i]) P[p++] = P[i];

    ConvexHull3D::getConvexHull(p);

    return 0;
}
```
