# 回文树/回文自动机（Palindrome Tree / Palindrome Automaton）

该实现的常数很大。

```c++
#include <cstdio>
#include <cstring>
#include <algorithm>

const int MAXN = 300005;
const int CHAR_SET = 26;

struct PalinT {
    struct Node {
        Node *c[CHAR_SET], *fail;
        int len, cnt;

        Node(int len = 0) : len(len), cnt(0), c(), fail(NULL) {}
    } *even, *odd, *last, _pool[MAXN], *_curr;
    char str[MAXN];
    int size;

    PalinT() : str() {
        _curr = _pool;
        last = even = new (_curr++) Node;
        even->fail = odd = new (_curr++) Node(-1);
        odd->fail = odd;
        str[size = 0] = -1;
    }

    Node *extend(int c) {
        str[++size] = c;
        Node *v = last;
        for (; str[size - v->len - 1] != str[size]; v = v->fail) {}

        Node *u = v->c[c];
        if (!u) {
            u = v->c[c] = new (_curr++) Node(v->len + 2);
            Node *p = v->fail;
            for (; str[size - p->len - 1] != str[size]; p = p->fail) {}
            u->fail = v == odd ? even : p->c[c];
        }
        u->cnt++;

        return last = u;
    }

    void build(char *begin, char *end) {
        for (char *p = begin; p != end; p++) extend(*p - 'a');
    }

    void count() {
        for (Node *p = _curr - 1; p >= _pool; p--) p->fail->cnt += p->cnt;
    }
} palinT;

int main() {
    static char s[MAXN];
    scanf("%s", s);

    int n = strlen(s);
    palinT.build(s, s + n);
    palinT.count();

    long long ans = 0;
    for (PalinT::Node *p = palinT._pool; p != palinT._curr; p++)
        ans = std::max(ans, (long long) p->len * p->cnt);

    printf("%lld\n", ans);

    return 0;
}
```
