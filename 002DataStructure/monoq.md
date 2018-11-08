# 单调队列（Monotonous Queue）

```c++
#include <cstdio>
#include <queue>

struct MonoQueue {
    std::deque<int> q, m;

    void push(int x) {
        q.push_back(x);
        while (!m.empty() && m.back() < x) m.pop_back();
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
} mq;

int main() {

    return 0;
}
```
