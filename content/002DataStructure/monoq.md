# 单调队列（Monotonous Queue）

```c++
#include <cstdio>
#include <queue>

template <typename T, typename Cmp = std::less<T> > // maximum by default
struct MonoQueue {
    std::deque<T> q, m;
    Cmp cmp;

    void push(int x) {
        q.push_back(x);
        while (!m.empty() && cmp(m.back(), x)) m.pop_back();
        m.push_back(x);
    }

    void pop() {
        int x = q.front();
        q.pop_front();
        if (x == m.front()) m.pop_front();
    }

    size_t size() {
        return q.size();
    }

    int top() {
        return m.front();
    }
};

int main() {

    return 0;
}
```
