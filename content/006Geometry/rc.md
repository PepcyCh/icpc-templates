# 旋转卡壳（Rotating Calipers ）

```c++
#include <cstdio>
#include <cfloat>
#include <cmath>
#include <algorithm>

const int MAXN = 50005;
const double EPS = 1e-9;

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

    friend double dist(const Point &a, const Point &b) {
        return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
    }

    void print() {
        if (std::abs(x) <= EPS) x = 0;
        if (std::abs(y) <= EPS) y = 0;
        printf("%.5f %.5f\n", x, y);
    }
    
    Point getPerpendicular() {
        double X = sqrt(1 / (1 + (x / y) * (x / y)));
        
        int sx, sy;
        if (x > 0 && y > 0) sx = 1, sy = -1;
        else if (x > 0 && y <= 0) sx = 1, sy = 1;
        else if (x <= 0 && y > 0) sx = -1, sy = -1;
        else sx = 1, sy = -1;
        
        return Point(sx * X, sy * std::sqrt(1 - X * X));
    }    
} P[MAXN], ch[MAXN];

struct Rectangle {
    Point p[4];
    
    void print() {
        int temp = 0;
        for (int i = 1; i < 4; i++)
            if (p[i].y < p[temp].y || (p[i].y == p[temp].y && p[i].x < p[temp].x)) temp = i;
        for (int i = 0; i < 4; i++) p[(i + temp) % 4].print();
    }
    
    Point &operator[](int i) { return p[i]; }
};

int getConvexHull(int n) {
    std::sort(P + 1, P + n + 1);
    
    int m = 0;
    for (int i = 1; i <= n; i++) {
        while (m > 1 && cross(ch[m - 1] - ch[m - 2], P[i] - ch[m - 2]) <= 0) m--;
        ch[m++] = P[i];
    }
    
    int k = m;
    for (int i = n; i; i--) {
        while (m > k && cross(ch[m - 1] - ch[m - 2], P[i] - ch[m - 2]) <= 0) m--;
        ch[m++] = P[i];
    }
    
    m > 1 ? m-- : 0;
    return m;
}

Rectangle ans;
double ansArea = DBL_MAX;
void rotatingCalipers(int n) {
    for (int curr = 0, right = 1, up = 1, left = 1; curr < n; curr++) {
        while (dot(ch[curr + 1] - ch[curr], ch[right + 1] - ch[right]) >= 0)
            right = (right + 1) % n;
        
        curr ? 0 : up = right;
        while (cross(ch[curr + 1] - ch[curr], ch[up + 1] - ch[up]) >= 0)
            up = (up + 1) % n;
        
        curr ? 0 : left = up;
        while (dot(ch[curr + 1] - ch[curr], ch[left + 1] - ch[left]) <= 0)
            left = (left + 1) % n;
        
        Point currV = ch[curr + 1] - ch[curr];
        double currLen = dist(ch[curr], ch[curr + 1]);
        double height = std::abs(cross(currV, ch[up] - ch[curr]) / currLen);
        double bottom = std::abs(dot(currV, ch[left] - ch[curr]) / currLen)
                      + std::abs(dot(currV, ch[right] - ch[curr]) / currLen);
                      
        double currArea = bottom * height;
        Point currPerpendicular = currV.getPerpendicular();
        
        if (currArea < ansArea) {
            ansArea = currArea;
            ans[0] = ch[curr] + currV * fabs((dot(currV, ch[right] - ch[curr])) / currLen) / currLen;
            ans[1] = ans[0] + currPerpendicular * height;
            ans[2] = ans[1] - currV * bottom / currLen;
            ans[3] = ans[2] - currPerpendicular * height;
        }
    }
}

int main() {
    int n;
    scanf("%d", &n);
    
    for (int i = 1; i <= n; i++) scanf("%lf %lf", &P[i].x, &P[i].y);
    
    int m = getConvexHull(n);
    rotatingCalipers(m);
    
    printf("%.5f\n", ansArea);
    ans.print();
    
    return 0;
}
```
