# 挑战多项式（多项式逆元、除法、平方根、ln、exp）

[LibreOJ #150](https://loj.ac/problem/150)

```c++
#include <cstdio>
#include <vector>
#include <tuple>
#include <algorithm>

const int MAXN = 262144 + 1;
const int MOD = 998244353;
const int G = 3;

long long qpow(long long a, long long n) {
    long long res = 1;
    for (; n; n >>= 1, a = a * a % MOD) if (n & 1) res = res * a % MOD;
    return res;
}

int numinv[MAXN];
void init() {
    numinv[1] = 1;
    for (int i = 2; i < MAXN; i++) numinv[i] = (long long) (MOD - MOD / i) * numinv[MOD % i] % MOD;
}

long long inv(long long x) {
    return x < MAXN ? numinv[x] : qpow(x, MOD - 2);
}

class NTT {
private:
    static const int N = 262144;

    long long omega[N + 1], omegaInv[N + 1];

    void init() {
        long long g = qpow(G, (MOD - 1) / N), ig = inv(g);
        omega[0] = omegaInv[0] = 1;
        for (int i = 1; i < N; i++) {
            omega[i] = omega[i - 1] * g % MOD;
            omegaInv[i] = omegaInv[i - 1] * ig % MOD;
        }
    }

    void reverse(long long *a, int n) const {
        for (int i = 0, j = 0; i < n; i++) {
            if (i < j) std::swap(a[i], a[j]);
            for (int l = n >> 1; (j ^= l) < l; l >>= 1) {}
        }
    }

    void transform(long long *a, int n, const long long *omega) const {
        reverse(a, n);

        for (int l = 2; l <= n; l <<= 1) {
            int hl = l >> 1;
            for (long long *x = a; x != a + n; x += l) {
                for (int i = 0; i < hl; i++) {
                    long long t = omega[N / l * i] * x[i + hl] % MOD;
                    x[i + hl] = (x[i] - t + MOD) % MOD;
                    x[i] += t;
                    x[i] >= MOD ? x[i] -= MOD : 0;
                }
            }
        }
    }

public:
    NTT() {
        init();
    }

    int extend(int n) const {
        int res = 1;
        while (res < n) res <<= 1;
        return res;
    }

    void dft(long long *a, int n) const {
        transform(a, n, omega);
    }

    void idft(long long *a, int n) const {
        transform(a, n, omegaInv);
        long long t = inv(n);
        for (int i = 0; i < n; i++) a[i] = a[i] * t % MOD;
    }
} ntt;

int modSqrt(int a, int p = MOD) {
    if (p == 2) return a % p;
    int x;
    if (qpow(a, (p - 1) >> 1) == 1) {
        if (p % 4 == 3) {
            x = qpow(a, (p + 1) >> 2);
        } else {
            long long w;
            for (w = 1; qpow((w * w - a + p) % p, (p - 1) >> 1) == 1; w++) {}
            long long b0 = w, b1 = 1;
            w = (w * w - a + p) % p;
            long long r0 = 1, r1 = 0;
            int exp = (p + 1) >> 1;
            for (; exp; std::tie(b0, b1) = std::make_tuple((b0 * b0 + b1 * b1 % p * w) % p, 2 * b0 * b1 % p), exp >>= 1) {
                if (exp & 1)
                    std::tie(r0, r1) = std::make_tuple((r0 * b0 + r1 * b1 % p * w) % p, (r0 * b1 + r1 * b0) % p);
            }
            x = r0;
        }
        if (x * 2 > p) x = p - x;
        return x;
    }
    return -1;
}

class Poly : public std::vector<long long> {
public:
    using std::vector<long long>::vector;

    static void dft(Poly &a, int n) {
        a.resize(n);
        ntt.dft(a.data(), n);
    }
    static void idft(Poly &a, int n) {
        a.resize(n);
        ntt.idft(a.data(), n);
    }

    Poly operator+(const Poly &rhs) const { // not commutative
        Poly res = *this;
        for (int i = 0; i < size(); i++) res[i] = (res[i] + rhs[i]) % MOD;
        return res;
    }

    Poly operator*(const Poly &rhs) const {
        if (size() < BRUTE_LIM || rhs.size() < BRUTE_LIM) {
            Poly res(size() + rhs.size() - 1);
            for (int i = 0; i < size(); i++)
                for (int j = 0; j < rhs.size(); j++)
                    res[i + j] = (res[i + j] + (*this)[i] * rhs[j]) % MOD;
            return res;
        }
        Poly t1 = *this, t2 = rhs;
        int n = t1.size() + t2.size() - 1;
        int N = ntt.extend(n);
        dft(t1, N);
        dft(t2, N);
        Poly res(N);
        for (int i = 0; i < N; i++) res[i] = t1[i] * t2[i] % MOD;
        idft(res, N);
        res.resize(n);
        return res;
    }

    static Poly inv(const Poly &a, int k = -1) {
        if (k == -1) k = a.size();
        if (k == 1) return { ::inv(a[0]) };
        Poly b = inv(a, (k + 1) >> 1), temp(a.begin(), a.begin() + k);
        int N = ntt.extend(2 * k - 1);
        dft(b, N);
        dft(temp, N);
        Poly res(N);
        for (int i = 0; i < N; i++) res[i] = (MOD + 2 - b[i] * temp[i] % MOD) * b[i] % MOD;
        idft(res, N);
        res.resize(k);
        return res;
    }

    static void div(const Poly &a, const Poly &b, Poly &d, Poly &r) {
        if (b.size() > a.size()) {
            d.clear();
            r = a;
            return;
        }

        int n = a.size(), m = b.size();

        Poly A = a, B = b;
        std::reverse(A.begin(), A.end());
        std::reverse(B.begin(), B.end());
        B.resize(n - m + 1);
        Poly iB = inv(B, n - m + 1);
        d = A * iB;
        d.resize(n - m + 1);
        std::reverse(d.begin(), d.end());

        r = b * d;
        r.resize(m - 1);
        for (int i = 0; i < m - 1; i++) r[i] = (a[i] - r[i] + MOD) % MOD;
    }

    static Poly derivative(const Poly &a) {
        Poly res(a.size() - 1);
        for (int i = 1; i < a.size(); i++) res[i - 1] = a[i] * i % MOD;
        return res;
    }

    static Poly integral(const Poly &a) {
        Poly res(a.size() + 1);
        for (int i = 0; i < a.size(); i++) res[i + 1] = a[i] * ::inv(i + 1) % MOD;
        return res;
    }

    static Poly log(const Poly &a) {
        Poly res = derivative(a) * inv(a);
        res.resize(a.size() - 1);
        return integral(res);
    }

    static Poly exp(const Poly &a, int k = -1) {
        if (k == -1) k = a.size();
        if (k == 1) return { 1 };
        Poly res = exp(a, (k + 1) >> 1);
        res.resize(k);
        Poly temp = log(res);
        for (auto &i : temp) i = i ? MOD - i : 0;
        ++temp[0];
        res = res * (temp + a);
        res.resize(k);
        return res;
    }

    static Poly sqrt(const Poly &a, int k = -1) {
        if (k == -1) k = a.size();
        if (k == 1) return { modSqrt(a[0]) };
        Poly res = sqrt(a, (k + 1) >> 1), temp(a.begin(), a.begin() + k);
        res.resize(k);
        res = res + temp * inv(res);
        res.resize(k);
        for (auto &i : res) i = (i % 2 ? ((i + MOD) >> 1) : (i >> 1));
        return res;
    }

    static Poly pow(const Poly &a, int n) {
        Poly res = log(a);
        for (auto &i : res) i = i * n % MOD;
        return exp(res);
    }


private:
    static const int BRUTE_LIM = 20;
};

int main() {
    init();

    int n, k;
    scanf("%d %d", &n, &k);

    Poly a(n + 1);
    for (int i = 0; i <= n; i++) scanf("%lld", &a[i]);

    return 0;
}
```
