# Sqrt Tree

在 O(n loglog n) 的时间内建树，在 O(1) 的时间内回答区间询问，在 O(sqrt(n)) 的时间内单点修改。

```c++
#include <cstdio>
#include <vector>
#include <algorithm>

inline int log2Up(int n) {
    int res = 0;
    while ((1 << res) < n) {
        res++;
    }
    return res;
}

template <typename T, T (*op)(T, T)>
class SqrtTree {
private:
    int n, lg, indexSz;
    std::vector<T> v;
    std::vector<int> clz, layers, onLayer;
    std::vector<std::vector<T> > pref, suf, between;

    void buildBlock(int layer, int l, int r) {
        pref[layer][l] = v[l];
        for (int i = l + 1; i < r; i++) {
            pref[layer][i] = op(pref[layer][i - 1], v[i]);
        }
        suf[layer][r - 1] = v[r - 1];
        for (int i = r - 2; i >= l; i--) {
            suf[layer][i] = op(v[i], suf[layer][i + 1]);
        }
    }

    void buildBetween(int layer, int lBound, int rBound, int betweenOffs) {
        int bSzLog = (layers[layer] + 1) >> 1;
        int bCntLog = layers[layer] >> 1;
        int bSz = 1 << bSzLog;
        int bCnt = (rBound - lBound + bSz - 1) >> bSzLog;
        for (int i = 0; i < bCnt; i++) {
            T ans;
            for (int j = i; j < bCnt; j++) {
                T add = suf[layer][lBound + (j << bSzLog)];
                ans = (i == j) ? add : op(ans, add);
                between[layer - 1][betweenOffs + lBound + (i << bCntLog) + j] = ans;
            }
        }
    }

    void buildBetweenZero() {
        int bSzLog = (lg + 1) >> 1;
        for (int i = 0; i < indexSz; i++) {
            v[n + i] = suf[0][i << bSzLog];
        }
        build(1, n, n + indexSz, (1 << lg) - n);
    }

    void updateBetweenZero(int bid) {
        int bSzLog = (lg + 1) >> 1;
        v[n + bid] = suf[0][bid << bSzLog];
        update(1, n, n + indexSz, (1 << lg) - n, n + bid);
    }

    void build(int layer, int lBound, int rBound, int betweenOffs) {
        if (layer >= (int)layers.size()) return;
        int bSz = 1 << ((layers[layer] + 1) >> 1);
        for (int l = lBound; l < rBound; l += bSz) {
            int r = std::min(l + bSz, rBound);
            buildBlock(layer, l, r);
            build(layer + 1, l, r, betweenOffs);
        }
        if (layer == 0) {
            buildBetweenZero();
        } else {
            buildBetween(layer, lBound, rBound, betweenOffs);
        }
    }

    void update(int layer, int lBound, int rBound, int betweenOffs, int x) {
        if (layer >= (int)layers.size()) return;
        int bSzLog = (layers[layer] + 1) >> 1;
        int bSz = 1 << bSzLog;
        int blockIdx = (x - lBound) >> bSzLog;
        int l = lBound + (blockIdx << bSzLog);
        int r = std::min(l + bSz, rBound);
        buildBlock(layer, l, r);
        if (layer == 0) {
            updateBetweenZero(blockIdx);
        } else {
            buildBetween(layer, lBound, rBound, betweenOffs);
        }
        update(layer + 1, l, r, betweenOffs, x);
    }

    T query(int l, int r, int betweenOffs, int base) {
        if (l == r) return v[l];
        if (l + 1 == r) return op(v[l], v[r]);
        int layer = onLayer[clz[(l - base) ^ (r - base)]];
        int bSzLog = (layers[layer] + 1) >> 1;
        int bCntLog = layers[layer] >> 1;
        int lBound = (((l - base) >> layers[layer]) << layers[layer]) + base;
        int lBlock = ((l - lBound) >> bSzLog) + 1;
        int rBlock = ((r - lBound) >> bSzLog) - 1;
        T ans = suf[layer][l];
        if (lBlock <= rBlock) {
            T add = (layer == 0)
                        ? (query(n + lBlock, n + rBlock, (1 << lg) - n, n))
                        : (between[layer - 1][betweenOffs + lBound + (lBlock << bCntLog) + rBlock]);
            ans = op(ans, add);
        }
        ans = op(ans, pref[layer][r]);
        return ans;
    }

public:
    inline T query(int l, int r) { return query(l, r, 0, 0); }

    inline void update(int x, const T &item) {
        v[x] = item;
        update(0, 0, n, 0, x);
    }

    SqrtTree(const std::vector<T> &a) : n((int) a.size()), lg(log2Up(n)), v(a), clz(1 << lg), onLayer(lg + 1) {
        clz[0] = 0;
        for (int i = 1; i < (int)clz.size(); i++) {
            clz[i] = clz[i >> 1] + 1;
        }
        int tlg = lg;
        while (tlg > 1) {
            onLayer[tlg] = (int)layers.size();
            layers.push_back(tlg);
            tlg = (tlg + 1) >> 1;
        }
        for (int i = lg - 1; i >= 0; i--) {
            onLayer[i] = std::max(onLayer[i], onLayer[i + 1]);
        }
        int betweenLayers = std::max(0, (int)layers.size() - 1);
        int bSzLog = (lg + 1) >> 1;
        int bSz = 1 << bSzLog;
        indexSz = (n + bSz - 1) >> bSzLog;
        v.resize(n + indexSz);
        pref.assign(layers.size(), std::vector<T>(n + indexSz));
        suf.assign(layers.size(), std::vector<T>(n + indexSz));
        between.assign(betweenLayers, std::vector<T>((1 << lg) + bSz));
        build(0, 0, n, 0);
    }
};

int Op(int a, int b) {
    return a + b; // gcd, min, max, &, |, ^
}

int main() {
    int n;
    scanf("%d", &n);
    std::vector<int> a(n);
    for (int i = 0; i < n; i++) scanf("%d", &a[i]);

    SqrtTree<int, Op> sqrtT(a);

    return 0;
}
```
