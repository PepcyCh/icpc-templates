# Pollard's Rho

```c++
long long gcd(long long a, long long b) {
    return b ? gcd(b, a % b) : a;
}

namespace PollardRho {
    long long g(long long x, long long n, long long c) {
        return (mul(x, x, n) + c) % n;
    }

    long long rho(long long n, long long c) {
        long long x = rand() % n, y = x, d = 1;

        for (long long i = 1, k = 2; d == 1; i++) {
            x = g(x, n, c);
            d = gcd(x > y ? x - y : y - x, n);
            if (x == y) return n;
            if (i == k) k <<= 1, y = x;
        }

        return d;
    }

    void find(long long n, long long c, std::map<long long, int> &res) {
        if (n == 1) return;
        if (isPrime(n)) {
            res[n]++;
            return;
        }

        long long p = n;
        while (p == n) p = rho(p, c++);
        find(p, c, res);
        find(n / p, c, res);
    }

    std::map<long long, int> divide(long long n) {
        static std::map<long long, int> res;
        res.clear();
        find(n, 1, res);
        return res;
    }
}
```
