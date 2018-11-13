# Miller-Rabin

```c++
long long pow(long long a, long long n, long long mod) {
    long long res = 1;
    for (; n; n >>= 1, a = mul(a, a, mod)) if (n & 1) res = mul(res, a, mod);
    return res;
}

bool isPrime(long long n) {
    const static int primes[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};

    long long s = 0, d = n - 1;
    while (d % 2 == 0) d /= 2, s++;
    if (s == 0) return n == 2;

    for (int i = 0; i < 12 && primes[i] < n; i++) {
        long long a = primes[i];
        if (pow(a, d, n) != 1) {
            bool flag = true;
            for (int r = 0; r < s; r++)
                if (flag && pow(a, d * (1 << r), n) == n - 1) flag = false;
            if (flag) return false;
        }
    }

    return true;
}
```
