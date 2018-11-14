# 三维偏序

```c++
#include <cstdio>
#include <algorithm>

const int MAXN = 100005;
const int MAXK = 200005;

struct Data {
    int a, b, c, cnt, ans;

    bool operator<(const Data &rhs) const {
        return a < rhs.a || (a == rhs.a && b < rhs.b) || (a == rhs.a && b == rhs.b && c < rhs.c);
    }

    bool operator==(const Data &rhs) const {
        return a == rhs.a && b == rhs.b && c == rhs.c;
    }
} a[MAXN];

struct BIT {
    int c[MAXK], n;

    static int lowbit(int x) {
        return x & -x;
    }

    void init(int n) {
        this->n = n;
    }

    void update(int pos, int d) {
        for (int i = pos; i <= n; i += lowbit(i)) c[i] += d;
    }

    int query(int pos) {
        int res = 0;
        for (int i = pos; i; i -= lowbit(i)) res += c[i];
        return res;
    }

    void clear(int pos) {
        for (int i = pos; i <= n; i += lowbit(i)) {
            if (c[i]) c[i] = 0;
            else break;
        }
    }
} bit;

void divide(Data *l, Data *r) {
    if (l == r) {
        l->ans += l->cnt - 1;
        return;
    }

    Data *mid = l + (r - l) / 2;
    divide(l, mid), divide(mid + 1, r);

    static Data temp[MAXN];
    for (Data *p = temp, *pl = l, *pr = mid + 1; p <= temp + (r - l); p++) {
        if (pr > r || (pl <= mid && pl->b <= pr->b)) {
            *p = *pl++;
            bit.update(p->c, p->cnt);
        } else {
            *p = *pr++;
            p->ans += bit.query(p->c);
        }
    }

    for (Data *p = temp, *q = l; q <= r; q++, p++) {
        bit.clear(q->c);
        *q = *p;
    }
}

int main() {
    int n, k;
    scanf("%d %d", &n, &k);

    for (int i = 0; i < n; i++) {
        scanf("%d %d %d", &a[i].a, &a[i].b, &a[i].c);
        a[i].cnt = 1;
    }

    std::sort(a, a + n);
    int cnt = 0;
    for (int i = 0; i < n; i++) {
        if (i && a[i] == a[i - 1]) a[cnt - 1].cnt++;
        else a[cnt++] = a[i];
    }

    bit.init(k);
    divide(a, a + cnt - 1);

    static int ans[MAXN];
    for (int i = 0; i < cnt; i++) ans[a[i].ans] += a[i].cnt;
    for (int i = 0; i < n; i++) printf("%d\n", ans[i]);

    return 0;
}
```
