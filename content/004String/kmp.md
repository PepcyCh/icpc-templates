# KMP

```c++
#include <cstdio>
#include <cstring>

const int MAXN = 1000005;

int fail[MAXN];
void calcFail(char *s) {
    int n = strlen(s + 1);

    fail[1] = 0;
    for (int i = 2; i <= n; i++) {
        int j = fail[i - 1];
        while (j && s[j + 1] != s[i]) j = fail[j];
        fail[i] = s[i] == s[j + 1] ? j + 1 : 0;
    }
}

int kmp(char *s, char *t) {
    int res = 0;
    int n = strlen(s + 1), m = strlen(t + 1);

    for (int i = 1, j = 0; i <= n; i++) {
        while (j && t[j + 1] != s[i]) j = fail[j];
        if (t[j + 1] == s[i]) j++;
        if (j == m) {
            res++;
            j = fail[j]; // j = 0 when not allowed overlapping
        }
    }

    return res;
}

int main() {
    static char a[MAXN], b[MAXN];
    scanf("%s %s", a + 1, b + 1);

    calcFail(b);

    printf("%d\n", kmp(a, b));

    return 0;
}
```
