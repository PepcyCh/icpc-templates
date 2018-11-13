# Manacher

```c++
#include <cstdio>
#include <cstring>
#include <algorithm>

const int MAXN = 1000005;

namespace Manacher {
    char s[MAXN << 1];
    int r[MAXN << 1], len;

    void prepare(char *str) {
        len = 0;
        s[++len] = '@';
        s[++len] = '#';
        int n = strlen(str);
        for (int i = 0; i < n; i++) s[++len] = str[i], s[++len] = '#';
        s[++len] = '\0';
    }

    void manacher() {
        int right = 0, pos = 0;
        for (int i = 1; i <= len; i++) {
            int x = right < i ? 1 : std::min(r[2 * pos - i], right - i);
            while (s[i + x] == s[i - x]) x++;
            r[i] = x;
            if (x + i > right) {
                right = x + i;
                pos = i;
            }
        }
    }

    void calc(char *str) {
        prepare(str);
        manacher();
    }
}

int main() {
    static char s[MAXN];
    scanf("%s", s);

    Manacher::calc(s);

    int ans = 0, len = Manacher::len, *r = Manacher::r;
    for (int i = 1; i <= len; i++) ans = std::max(ans, r[i] - 1);

    printf("%d\n", ans);

    return 0;
}
```
