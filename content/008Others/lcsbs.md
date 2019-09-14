# bitset 优化最长公共子序列

```c++
#include <cstdio>
#include <cstring>
#include <algorithm>

const int MAXN = 5005, SIGMA = 26;

struct Bitset {
    static constexpr int W = 62;

    static constexpr size_t size = MAXN / W + !!(MAXN % W);
    unsigned long long u[size];

    Bitset &operator=(const Bitset &rhs) {
        std::copy(rhs.u, rhs.u + size, u);
        return *this;
    }

    void reset() {
        std::fill(u, u + size, 0);
    }
    void set(int x) {
        u[x / W] |= 1ull << (x % W);
    }
    bool get(int x) const {
        return u[x / W] & (1ull << (x % W));
    }

    Bitset operator|(const Bitset &rhs) const {
        Bitset res;
        for (int i = 0; i < size; i++) res.u[i] = u[i] | rhs.u[i];
        return res;
    }
    Bitset operator&(const Bitset &rhs) const {
        Bitset res;
        for (int i = 0; i < size; i++) res.u[i] = u[i] & rhs.u[i];
        return res;
    }
    Bitset operator^(const Bitset &rhs) const {
        Bitset res;
        for (int i = 0; i < size; i++) res.u[i] = u[i] ^ rhs.u[i];
        return res;
    }

    Bitset operator-(const Bitset &rhs) const {
        Bitset res;
        for (int i = 0; i < size; i++) res.u[i] = u[i] - rhs.u[i];
        for (int i = 0; i < size; ++ i) if (res.u[i] < 0) {
            res.u[i] += 1ll << W;
            --res.u[i + 1];
        }
        res.u[size] = 0;
        return res;
    }

    void shlOr1() {
        unsigned long long c = 1;
        for (int i = 0; i < size; ++ i) {
            u[i] <<= 1; u[i] |= c;
            c = u[i] >> W & 1;
            u[i] ^= c << W;
        }
    }

    size_t count() const {
        size_t res = 0;
        for (int i = 0; i < size; ++ i) res += __builtin_popcountll(u[i]);
        return res;
    }
};
Bitset row, as[SIGMA], x;

char s[MAXN], t[MAXN];

int main() {
    scanf("%s %s", s, t);
    int n = strlen(s), m = strlen(t);

    for (int i = 0; i < SIGMA; i++) as[i].reset();
    for (int i = 0; i < n; i++) as[s[i] - 'a'].set(i);

    row.reset();
    for (int i = 0; i < m; i++) {
        x = row | as[t[i] - 'a'];
        row.shlOr1();
        row = x - row;
        row = (row ^ x) & x;
        printf("%d\n", row.count());
    }
    return 0;
}
```
