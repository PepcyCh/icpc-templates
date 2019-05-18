# 多项式逆元

若不需要写拆系数 FFT，或可以使用 NTT，则两次多项式乘法可用一次代替。

```c++
#include <cstdio>
#include <cmath>
#include <algorithm>

const int MAXN = 262144 + 1;
const double PI = std::acos(-1.0);
const int MOD = 1000000007;

struct Complex {
    double r, i;

    Complex(double r = 0, double i = 0) : r(r), i(i) {}

    Complex conj() const { return Complex(r, -i); }

    Complex operator-(const Complex &rhs) const { return Complex(r - rhs.r, i - rhs.i); }
    Complex operator+(const Complex &rhs) const { return Complex(r + rhs.r, i + rhs.i); }
    Complex operator*(const Complex &rhs) const { return Complex(r * rhs.r - i * rhs.i, r * rhs.i + i * rhs.r); }
    Complex operator/(double rhs) const { return Complex(r / rhs, i / rhs); }
};

namespace FFT {
    const int N = 262144;

    Complex omega[::MAXN], omegaInv[::MAXN];

    void init() {
        double per = 2 * PI / N;
        for (int i = 0; i < N; i++) {
            omega[i] = Complex(std::cos(i * per), std::sin(i * per));
            omegaInv[i] = omega[i].conj();
        }
    }

    int extend(int n) {
        int res = 1;
        while (res < n) res <<= 1;
        return res;
    }

    void reverse(Complex *a, int n) {
        for (int i = 0, j = 0; i < n; i++) {
            if (i < j) std::swap(a[i], a[j]);
            for (int l = n >> 1; (j ^= l) < l; l >>= 1) {}
        }
    }

    void transform(Complex *a, int n, Complex *omega) {
        reverse(a, n);

        for (int l = 2; l <= n; l <<= 1) {
            int hl = l >> 1;
            for (Complex *x = a; x != a + n; x += l) {
                for (int i = 0; i < hl; i++) {
                    Complex t = omega[N / l * i] * x[i + hl];
                    x[i + hl] = x[i] - t;
                    x[i] = x[i] + t;
                }
            }
        }
    }

    void dft(Complex *a, int n) {
        transform(a, n, omega);
    }

    void idft(Complex *a, int n) {
        transform(a, n, omegaInv);
        for (int i = 0; i < n; i++) a[i] = a[i] / n;
    }
}

void modMul(long long *a, long long *b, int n, long long *res) {
    static Complex a0[MAXN], a1[MAXN], b0[MAXN], b1[MAXN];
    static const int M = (1 << 15) - 1;

    for (int i = 0; i < n; i++) {
        a0[i] = a[i] >> 15;
        a1[i] = a[i] & M;
        b0[i] = b[i] >> 15;
        b1[i] = b[i] & M;
    }
    FFT::dft(a0, n), FFT::dft(a1, n);
    FFT::dft(b0, n), FFT::dft(b1, n);
    for (int i = 0; i < n; i++) {
        Complex _a = a0[i], _b = a1[i], _c = b0[i], _d = b1[i];
        a0[i] = _a * _c;
        a1[i] = _a * _d + _b * _c;
        b0[i] = _b * _d;
    }
    FFT::idft(a0, n), FFT::idft(a1, n), FFT::idft(b0, n);
    for (int i = 0; i < n; i++) {
        res[i] = ((((long long) (a0[i].r + 0.5) % MOD) << 30) % MOD
                + (((long long) (a1[i].r + 0.5) % MOD) << 15) % MOD
                  + (long long) (b0[i].r + 0.5) % MOD) % MOD;
    }
}

void polyInverse(long long *a, long long *res, int k) {
    if (k == 1) {
        res[0] = 1;
        return;
    }
    polyInverse(a, res, (k + 1) >> 1);

    static long long t1[MAXN], t2[MAXN];
    int N = FFT::extend(k << 1);
    std::copy(a, a + k, t1);
    std::fill(t1 + k, t1 + N, 0);
    modMul(res, res, N, t2);
    modMul(t1, t2, N, t1);
    for (int i = 0; i < k; i++) res[i] = (2 * res[i] % MOD - t1[i] + MOD) % MOD;
    std::fill(res + k, res + N, 0);
}

long long a[MAXN], res[MAXN];

int main() {
    FFT::init();

    int n, k;
    scanf("%d %d", &n, &k);
    for (int i = 0; i <= n; i++) scanf("%lld", &a[i]);

    polyInverse(a, res, k);
    for (int i = 0; i < k; i++) printf("%lld%c", res[i], " \n"[i == k - 1]);

    return 0;
}
```
