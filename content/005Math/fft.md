# 快速傅立叶变换（Fast Fourier Transformation）

### FFT

```c++
#include <cstdio>
#include <cmath>
#include <algorithm>

const int MAXN = 100005;
const int MAXN_EXTEND = 262144;
const double PI = std::acos(-1.0);

struct Complex {
    double r, i;

    Complex(double r = 0, double i = 0) : r(r), i(i) {}

    Complex conj() const {
        return Complex(r, -i);
    }

    Complex operator+(const Complex &rhs) const {
        return Complex(r + rhs.r, i + rhs.i);
    }
    Complex &operator+=(const Complex &rhs) {
        return *this = *this + rhs;
    }

    Complex operator-(const Complex &rhs) const {
        return Complex(r - rhs.r, i - rhs.i);
    }

    Complex operator*(const Complex &rhs) const {
        return Complex(r * rhs.r - i * rhs.i, r * rhs.i + i * rhs.r);
    }
    Complex &operator*=(const Complex &rhs) {
        return *this = *this * rhs;
    }

    Complex operator/(const double rhs) const {
        return Complex(r / rhs, i / rhs);
    }
    Complex &operator/=(const double rhs) {
        return *this = *this / rhs;
    }
};

namespace FFT {
    static const int N = 262144;

    Complex omega[N], omegaInv[N];

    void init() {
        for (int i = 0; i < N; i++) {
            omega[i] = Complex(std::cos(2 * PI * i / N), std::sin(2 * PI * i / N));
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
                    x[i] += t;
                }
            }
        }
    }

    void dft(Complex *a, int n) {
        transform(a, n, omega);
    }

    void idft(Complex *a, int n) {
        transform(a, n, omegaInv);
        for (int i = 0; i < n; i++) a[i] /= n;
    }
}

int main() {
    int n, m;
    scanf("%d %d", &n, &m);

    static Complex a[MAXN_EXTEND], b[MAXN_EXTEND];
    for (int i = 0; i <= n; i++) scanf("%lf", &a[i].r);
    for (int i = 0; i <= m; i++) scanf("%lf", &b[i].r);

    FFT::init();
    int N = FFT::extend(n + m + 1);

    FFT::dft(a, N);
    FFT::dft(b, N);
    for (int i = 0; i < N; i++) a[i] *= b[i];
    FFT::idft(a, N);

    for (int i = 0; i < n + m + 1; i++)
        printf("%.0lf%c", a[i].r + 0.001, " \n"[i == n + m]);

    return 0;
}
```

### 两次 DFT 的多项式乘法

```c++
#include <cstdio>
#include <cmath>
#include <algorithm>

const int MAXN = 262144;
const double PI = std::acos(-1.0);

struct Complex {
    double r, i;

    Complex(double r = 0, double i = 0) : r(r), i(i) {}

    Complex conj() const {
        return Complex(r, -i);
    }

    Complex operator+(const Complex &rhs) const {
        return Complex(r + rhs.r, i + rhs.i);
    }
    Complex &operator+=(const Complex &rhs) {
        return *this = *this + rhs;
    }

    Complex operator-(const Complex &rhs) const {
        return Complex(r - rhs.r, i - rhs.i);
    }

    Complex operator*(const Complex &rhs) const {
        return Complex(r * rhs.r - i * rhs.i, r * rhs.i + i * rhs.r);
    }
    Complex &operator*=(const Complex &rhs) {
        return *this = *this * rhs;
    }

    Complex operator/(const double rhs) const {
        return Complex(r / rhs, i / rhs);
    }
    Complex &operator/=(const double rhs) {
        return *this = *this / rhs;
    }
};

namespace FFT {
    static const int N = 262144;

    Complex omega[N], omegaInv[N];

    void init() {
        for (int i = 0; i < N; i++) {
            omega[i] = Complex(std::cos(2 * PI * i / N), std::sin(2 * PI * i / N));
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
                    x[i] += t;
                }
            }
        }
    }

    void dft(Complex *a, int n) {
        transform(a, n, omega);
    }

    void idft(Complex *a, int n) {
        transform(a, n, omegaInv);
        for (int i = 0; i < n; i++) a[i] /= n;
    }
}

int main() {
    FFT::init();

    int n, m;
    scanf("%d %d", &n, &m);

    static Complex a[MAXN], b[MAXN];

    for (int i = 0, x; i <= n; i++) {
        scanf("%d", &x);
        a[i] = Complex(x, 0);
    }

    for (int i = 0, x; i <= m; i++) {
        scanf("%d", &x);
        a[i] += Complex(0, x);
    }

    int len = n + m + 1;
    int N = FFT::extend(len);

    FFT::dft(a, N);
    for (int i = 1; i < N; i++) {
        double x1 = a[i].r, y1 = a[i].i;
        double x2 = a[N - i].r, y2 = a[N - i].i;
        Complex t1((x1 + x2) * 0.5, (y1 - y2) * 0.5);
        Complex t2((y1 + y2) * 0.5, (x2 - x1) * 0.5);
        b[i] = t1 * t2;
    }
    b[0] = a[0].r * a[0].i;
    FFT::idft(b, N);

    for (int i = 0; i < len; i++)
        printf("%d%c", (int) (round(b[i].r)), " \n"[i == len - 1]);

    return 0;
}
```

### 任意模数 FFT

```c++
#include <cstdio>
#include <cmath>
#include <algorithm>

const int MAXN = 262144;
const int MOD = 1000000007;
const long double PI = std::acos(-1.0);

struct Complex {
    long double r, i;

    Complex(long double r = 0, long double i = 0) : r(r), i(i) {}

    Complex conj() const {
        return Complex(r, -i);
    }

    Complex operator+(const Complex &rhs) const {
        return Complex(r + rhs.r, i + rhs.i);
    }
    Complex &operator+=(const Complex &rhs) {
        return *this = *this + rhs;
    }

    Complex operator-(const Complex &rhs) const {
        return Complex(r - rhs.r, i - rhs.i);
    }
    Complex &operator-=(const Complex &rhs) {
        return *this = *this - rhs;
    }

    Complex operator*(const Complex &rhs) const {
        return Complex(r * rhs.r - i * rhs.i, r * rhs.i + i * rhs.r);
    }
    Complex &operator*=(const Complex &rhs) {
        return *this = *this * rhs;
    }

    Complex operator/(const long double rhs) const {
        return Complex(r / rhs, i / rhs);
    }
    Complex &operator/=(const long double rhs) {
        return *this = *this / rhs;
    }
};

namespace FFT {
    static const int N = 262144;

    Complex omega[N], omegaInv[N];

    void init() {
        for (int i = 0; i < N; i++) {
            omega[i] = Complex(std::cos(2 * PI * i / N), std::sin(2 * PI * i / N));
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
                    x[i] += t;
                }
            }
        }
    }

    void dft(Complex *a, int n) {
        transform(a, n, omega);
    }

    void idft(Complex *a, int n) {
        transform(a, n, omegaInv);
        for (int i = 0; i < n; i++) a[i] /= n;
    }
}

void lowAccuracyModMul(int *a, int *b, int n, int *res) {
    static Complex A[MAXN], B[MAXN];
    for (int i = 0; i < n; i++) A[i] = Complex(a[i], b[i]);
    FFT::dft(A, n);
    for (int i = 1; i < n; i++) {
        long double x1 = A[i].r, y1 = A[i].i;
        long double x2 = A[n - i].r, y2 = A[n - i].i;
        Complex t1((x1 + x2) * 0.5, (y1 - y2) * 0.5);
        Complex t2((y1 + y2) * 0.5, (x2 - x1) * 0.5);
        B[i] = t1 * t2;
    }
    B[0] = A[0].r * A[0].i;
    FFT::idft(B, n);
    for (int i = 0; i < n; i++)
        res[i] = ((long long) (B[i].r + 0.5)) % MOD;
}

void modMul(int *a, int *b, int n, int *res) {
    static int a0[MAXN], a1[MAXN], b0[MAXN], b1[MAXN];
    static const int M = 1000;

    for (int i = 0; i < n; i++) {
        a0[i] = a[i] / M;
        b0[i] = b[i] / M;
    }
    lowAccuracyModMul(a0, b0, n, a0);
    for (int i = 0; i < n; i++) {
        res[i] = 1ll * a0[i] * M * M % MOD;
        a1[i] = a[i] % M;
        b1[i] = b[i] % M;
    }
    lowAccuracyModMul(a1, b1, n, a1);
    for (int i = 0; i < n; i++) {
        res[i] += a1[i];
        res[i] >= MOD ? res[i] -= MOD : 0;

        a0[i] = (a0[i] + a1[i]) % MOD;
        a1[i] = a[i] / M + a[i] % M;
        b1[i] = b[i] / M + b[i] % M;
    }
    lowAccuracyModMul(a1, b1, n, a1);
    for (int i = 0; i < n; i++)
        res[i] = (1ll * M * (a1[i] - a0[i] + MOD) + res[i]) % MOD;
}

int main() {
    FFT::init();

    int n, m;
    scanf("%d %d", &n, &m);

    static int a[MAXN], b[MAXN], ans[MAXN];
    for (int i = 0; i <= n; i++) scanf("%d", &a[i]);
    for (int i = 0; i <= m; i++) scanf("%d", &b[i]);
    for (int i = 0; i <= n; i++) a[i] < 0 ? a[i] += MOD : 0;
    for (int i = 0; i <= m; i++) b[i] < 0 ? b[i] += MOD : 0;

    int len = n + m + 1;
    int p = FFT::extend(len);
    modMul(a, b, p, ans);
    for (int i = 0; i < len; i++) printf("%d%c", ans[i], " \n"[i == len - 1]);

    return 0;
}
```

### NTT (Number Theoretic Transformation)

```c++
#include <cstdio>
#include <algorithm>

const int MAXN = 100005;
const int MAXN_EXTEND = 262144;
const int MOD = 998244353;
const int G = 3;

long long pow(long long a, long long n) {
    long long res = 1;
    for (; n; n >>= 1, a = a * a % MOD) if (n & 1) res = res * a % MOD;
    return res;
}

void exgcd(long long a, long long b, long long &x, long long &y) {
    if (!b) x = 1, y = 0;
    else exgcd(b, a % b, y, x), y -= x * (a / b);
}

long long inv(long long x) {
    long long res, temp;
    exgcd(x, MOD, res, temp);
    return (res % MOD + MOD) % MOD;
}

namespace NTT {
    static const int N = 262144;

    long long omega[N], omegaInv[N];

    void init() {
        long long g = pow(G, (MOD - 1) / N), ig = inv(g);
        omega[0] = omegaInv[0] = 1;
        for (int i = 1; i < N; i++) {
            omega[i] = omega[i - 1] * g % MOD;
            omegaInv[i] = omegaInv[i - 1] * ig % MOD;
        }
    }

    int extend(int n) {
        int res = 1;
        while (res < n) res <<= 1;
        return res;
    }

    void reverse(long long *a, int n) {
        for (int i = 0, j = 0; i < n; i++) {
            if (i < j) std::swap(a[i], a[j]);
            for (int l = n >> 1; (j ^= l) < l; l >>= 1) {}
        }
    }

    void transform(long long *a, int n, long long *omega) {
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

    void dft(long long *a, int n) {
        transform(a, n, omega);
    }

    void idft(long long *a, int n) {
        transform(a, n, omegaInv);
        long long t = inv(n);
        for (int i = 0; i < n; i++) (a[i] *= t) %= MOD;
    }
}

int main() {
    int n, m;
    scanf("%d %d", &n, &m);

    static long long a[MAXN_EXTEND], b[MAXN_EXTEND];
    for (int i = 0; i <= n; i++) scanf("%lld", &a[i]);
    for (int i = 0; i <= m; i++) scanf("%lld", &b[i]);

    NTT::init();
    int N = NTT::extend(n + m + 1);

    NTT::dft(a, N);
    NTT::dft(b, N);
    for (int i = 0; i < N; i++) (a[i] *= b[i]) %= MOD;
    NTT::idft(a, N);

    for (int i = 0; i < n + m + 1; i++)  printf("%lld%c", a[i], " \n"[i == n + m]);

    return 0;
}
```

### NTT 模数及原根表

| $$r \dot 2 ^ k + 1$$ | $$r$$ | $$k$$ | $$g$$ |
|:---|:---|:---|:---|
| $$3$$ | $$1$$ | $$1$$ | $$2$$ | 
| $$5$$ | $$1$$ | $$2$$ | $$2$$ | 
| $$17$$ | $$1$$ | $$4$$ | $$3$$ | 
| $$97$$ | $$3$$ | $$5$$ | $$5$$ | 
| $$193$$ | $$3$$ | $$6$$ | $$5$$ | 
| $$257$$ | $$1$$ | $$8$$ | $$3$$ | 
| $$7681$$ | $$15$$ | $$9$$ | $$17$$ | 
| $$12289$$ | $$3$$ | $$12$$ | $$11$$ | 
| $$40961$$ | $$5$$ | $$13$$ | $$3$$ | 
| $$65537$$ | $$1$$ | $$16$$ | $$3$$ | 
| $$786433$$ | $$3$$ | $$18$$ | $$10$$ | 
| $$5767169$$ | $$11$$ | $$19$$ | $$3$$ | 
| $$7340033$$ | $$7$$ | $$20$$ | $$3$$ | 
| $$23068673$$ | $$11$$ | $$21$$ | $$3$$ | 
| $$104857601$$ | $$25$$ | $$22$$ | $$3$$ | 
| $$167772161$$ | $$5$$ | $$25$$ | $$3$$ | 
| $$469762049$$ | $$7$$ | $$26$$ | $$3$$ | 
| $$998244353$$ | $$110$$ | $$23$$ | $$3$$ | 
| $$1004535809$$ | $$479$$ | $$21$$ | $$3$$ | 
| $$2013265921$$ | $$15$$ | $$27$$ | $$31$$ | 
| $$2281701377$$ | $$17$$ | $$27$$ | $$3$$ | 
| $$3221225473$$ | $$3$$ | $$30$$ | $$5$$ | 
| $$75161927681$$ | $$35$$ | $$31$$ | $$3$$ | 
| $$77309411329$$ | $$9$$ | $$33$$ | $$7$$ | 
| $$206158430209$$ | $$3$$ | $$36$$ | $$22$$ | 
| $$2061584302081$$ | $$15$$ | $$37$$ | $$7$$ | 
