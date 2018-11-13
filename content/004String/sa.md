# 后缀数组（Suffix Array）

```c++
#include <cstdio>
#include <cstring>
#include <algorithm>

const int MAXN = 1000005;

namespace SuffixArray {
    int n, sa[MAXN], rank[MAXN], height[MAXN];
    char str[MAXN];

    void buildSA(int m) {
        static int fir[MAXN], sec[MAXN], buc[MAXN], temp[MAXN];
        n = strlen(str);
        str[n++] = 0;

        std::fill(buc, buc + m, 0);
        for (int i = 0; i < n; i++) buc[(int) str[i]]++;
        for (int i = 1; i < m; i++) buc[i] += buc[i - 1];
        for (int i = 0; i < n; i++) rank[i] = buc[(int) str[i]] - 1;

        for (int l = 1; l < n; l <<= 1) {
            for (int i = 0; i < n; i++)
                fir[i] = rank[i], sec[i] = i + l < n ? rank[i + l] : 0;

            std::fill(buc, buc + n, 0);
            for (int i = 0; i < n; i++) buc[sec[i]]++;
            for (int i = 1; i < n; i++) buc[i] += buc[i - 1];
            for (int i = n - 1; ~i; i--) temp[--buc[sec[i]]] = i;

            std::fill(buc, buc + n, 0);
            for (int i = 0; i < n; i++) buc[fir[i]]++;
            for (int i = 1; i < n; i++) buc[i] += buc[i - 1];
            for (int i = n - 1; ~i; i--) sa[--buc[fir[temp[i]]]] = temp[i];

            rank[sa[0]] = 0;
            bool unique = true;
            for (int i = 1; i < n; i++) {
                rank[sa[i]] = rank[sa[i - 1]];
                if (fir[sa[i]] == fir[sa[i - 1]] && sec[sa[i]] == sec[sa[i - 1]])
                    unique = false;
                else rank[sa[i]]++;
            }

            if (unique) break;
        }
    }

    void getHeight() {
        int k = 0;
        for (int i = 0; i < n - 1; i++) {
            k ? k-- : 0;
            int j = sa[rank[i] - 1];
            while (str[i + k] == str[j + k]) k++;
            height[rank[i]] = k;
        }
    }
}

int main() {
    char *str = SuffixArray::str;
    scanf("%s", str);

    SuffixArray::buildSA(128);
    int *sa = SuffixArray::sa, n = SuffixArray::n;
    for (int i = 1; i < n; i++) printf("%d%c", sa[i] + 1, " \n"[i == n - 1]);

    SuffixArray::getHeight();
    int *height = SuffixArray::height + 1;
    for (int i = 1; i < n - 1; i++) printf("%d%c", height[i], " \n"[i == n - 2]);

    return 0;
}
```
