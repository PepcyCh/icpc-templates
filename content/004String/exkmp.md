# 扩展 KMP

```c++
#include <cstdio>
#include <cstring>

const int MAXN = 1000005;

namespace ExKmp {
    int next[MAXN], extend[MAXN];

    void getNext(char *str, int n) {
        int i = 0;
        while (str[i] == str[i + 1] && i + 1 < n) i++;
        next[0] = n;
        next[1] = i;
        int pos = 1;
        for (i = 2; i < n; i++) {
            if (next[i - pos] + i < next[pos] + pos) next[i] = next[i - pos];
            else {
                int j = next[pos] + pos - i;
                if (j < 0) j = 0;
                while (i + j < n && str[j] == str[j + i]) ++j;
                next[i] = j;
                pos = i;
            }
        }
    }

    void getExtand(char *s, char *t) {
        int n = strlen(s), m = strlen(t);
        getNext(t, m);
        int i = 0;
        while (s[i] == t[i] && i < m && i < n) i++;
        extend[0] = i;
        int pos = 0;
        for (i = 1; i < n; i++) {
            if (next[i - pos] + i < extend[pos] + pos) extend[i] = next[i - pos];
            else {
                int j = extend[pos] + pos - i;
                if (j < 0) j = 0;
                while (i + j < n && j < m && s[j + i] == t[j]) ++j;
                extend[i] = j;
                pos = i;
            }
        }
    }
}

char s[MAXN], t[MAXN];

int main() {
    scanf("%s %s", s, t);
    ExKmp::getExtand(s, t);

    return 0;
}
```
