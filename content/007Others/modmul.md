# $$O(1)$$ 取模乘

```c++
long long mul(long long a, long long b, long long p) {
    return (a * b - (long long) (a / (long double) p * b + 1e-3) * p + p) % p;
}
```
