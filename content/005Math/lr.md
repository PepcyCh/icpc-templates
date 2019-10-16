# 常系数线性齐次递推

* 暴力多项式取模
* FFT 多项式取模

### 暴力多项式取模

```c++
#include <cstdio>
#include <algorithm>

const int MOD = 1000000007;
const int MAXK = 2005;

long long a[MAXK];
void modMul(long long *A, long long *B, int n, long long *r) {
    static long long res[MAXK << 1];
    std::fill(res, res + (n << 1), 0);

    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
        res[i + j] += A[i] * B[j] % MOD;
        res[i + j] >= MOD ? res[i + j] -= MOD : 0;
    }

    for (int i = (n << 1) - 2; i >= n; i--) if (res[i]) for (int j = n - 1; ~j; j--) {
        res[i - n + j] += res[i] * a[j] % MOD;
        res[i - n + j] >= MOD ? res[i - n + j] -= MOD : 0;
    }

    for (int i = 0; i < n; i++) r[i] = res[i];
}

int main() {
    int n, k;
    scanf("%d %d", &n, &k);

    static long long h[MAXK];
    for (int i = 0; i < k; i++) {
        scanf("%lld", &a[i]);
        a[i] < 0 ? a[i] += MOD : 0;
    }
    std::reverse(a, a + k);
    a[k] = 1;

    for (int i = 0; i < k; i++) {
        scanf("%lld", &h[i]);
        h[i] < 0 ? h[i] += MOD : 0;
    }

    if (n < k) return printf("%lld\n", h[n]), 0;

    static long long m[MAXK], t[MAXK];
    m[0] = t[1] = 1;

    for (int i = n; i; i >>= 1, modMul(t, t, k, t)) if (i & 1) modMul(m, t, k, m);

    long long hn = 0;
    for (int i = 0; i < k; i++) hn = (hn + m[i] * h[i] % MOD) % MOD;

    printf("%lld\n", hn);

    return 0;
}
```

### FFT 多项式取模

```c++
#include <cstdio>
#include <cmath>
#include <algorithm>

const int MAXN = 131072 + 1;
const int MOD = 1000000007;
const double PI = std::acos(-1);

long long qpow(long long a, long long n) {
    long long res = 1;
    for (; n; n >>= 1, a = a * a % MOD) if (n & 1) res = res * a % MOD;
    return res;
}

long long inv(long long x) {
    return qpow(x, MOD - 2);
}

struct Complex {
    double r, i;

    Complex(double r = 0, double i = 0) : r(r), i(i) {}

    Complex conj() const { return Complex(r, -i); }

    Complex operator-(const Complex &rhs) const { return Complex(r - rhs.r, i - rhs.i); }
    Complex operator+(const Complex &rhs) const { return Complex(r + rhs.r, i + rhs.i); }
    Complex operator*(const Complex &rhs) const { return Complex(r * rhs.r - i * rhs.i, r * rhs.i + i * rhs.r); }
    Complex operator/(double rhs) const { return Complex(r / rhs, i / rhs); }
};

class FFT {
private:
    static const int N = 131072;

    Complex omega[N + 1], omegaInv[N + 1];

    void init() {
        double per = 2 * PI / N;
        for (int i = 0; i < N; i++) {
            omega[i] = Complex(std::cos(i * per), std::sin(i * per));
            omegaInv[i] = omega[i].conj();
        }
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

public:
    FFT() { init(); }

    int extend(int n) {
        int res = 1;
        while (res < n) res <<= 1;
        return res;
    }

    void dft(Complex *a, int n) {
        transform(a, n, omega);
    }

    void idft(Complex *a, int n) {
        transform(a, n, omegaInv);
        for (int i = 0; i < n; i++) a[i] = a[i] / n;
    }
} fft;

void polyMul(int n, int m, long long *a, long long *b, long long *res) {
    static Complex a0[MAXN], a1[MAXN], b0[MAXN], b1[MAXN];
    static const int M = (1 << 15) - 1;
    int N = fft.extend(n + m - 1);

    for (int i = 0; i < n; i++) {
        a0[i] = a[i] >> 15;
        a1[i] = a[i] & M;
    }
    std::fill(a0 + n, a0 + N, 0);
    std::fill(a1 + n, a1 + N, 0);
    for (int i = 0; i < m; i++) {
        b0[i] = b[i] >> 15;
        b1[i] = b[i] & M;
    }
    std::fill(b0 + n, b0 + N, 0);
    std::fill(b1 + n, b1 + N, 0);

    fft.dft(a0, N), fft.dft(a1, N);
    fft.dft(b0, N), fft.dft(b1, N);
    for (int i = 0; i < N; i++) {
        Complex _a = a0[i], _b = a1[i], _c = b0[i], _d = b1[i];
        a0[i] = _a * _c;
        a1[i] = _a * _d + _b * _c;
        b0[i] = _b * _d;
    }
    fft.idft(a0, N), fft.idft(a1, N), fft.idft(b0, N);
    for (int i = 0; i < N; i++) {
        res[i] = ((((long long) (a0[i].r + 0.5) % MOD) << 30) % MOD
                + (((long long) (a1[i].r + 0.5) % MOD) << 15) % MOD
                  + (long long) (b0[i].r + 0.5) % MOD) % MOD;
    }
}

void polyInv(int k, long long *a, long long *res) {
    if (k == 1) {
        res[0] = inv(a[0]);
        return;
    }
    polyInv((k + 1) >> 1, a, res);

    static long long t1[MAXN], t2[MAXN];
    int N = fft.extend(k << 1);
    std::copy(a, a + k, t1);
    std::fill(t1 + k, t1 + N, 0);
    polyMul(N, N, res, res, t2);
    polyMul(N, N, t1, t2, t1);
    for (int i = 0; i < k; i++) res[i] = (2 * res[i] % MOD - t1[i] + MOD) % MOD;
    std::fill(res + k, res + N, 0);
}

void polyDiv(int n, int m, long long *A, long long *B, long long *D, long long *R) {
	static long long A0[MAXN], B0[MAXN];

	int t = n - m + 1;
    int N = fft.extend(t << 1);

    std::fill(A0, A0 + N, 0);
    std::reverse_copy(B, B + m, A0);
	polyInv(t, A0, B0);

    std::reverse_copy(A, A + n, A0);

    polyMul(t, t, A0, B0, A0);

    std::reverse(A0, A0 + t);
    std::copy(A0, A0 + t, D);

    N = fft.extend(n);
    std::copy(B, B + m, B0);
    polyMul(t, m, A0, B0, A0);
	for (int i = 0; i < m; i++) R[i] = (A[i] - A0[i] + MOD) % MOD;
    std::fill(R + m, R + N, 0);
}

void polyPow(int n, long long *a, long long *mod, int k, long long *res) {
    res[0] = 1;
    static long long D[MAXN], R[MAXN];
    for (; k; k >>= 1) {
        if (k & 1) {
            polyMul(n, n, res, a, res);
            polyDiv(2 * n - 1, n, res, mod, D, R);
            std::copy(R, R + n, res);
        }
        polyMul(n, n, a, a, a);
        polyDiv(2 * n - 1, n, a, mod, D, R);
        std::copy(R, R + n, a);
    }
}

long long a[MAXN], x[MAXN], f[MAXN], t[MAXN], m[MAXN];

int main() {
    int n, k;
    scanf("%d %d", &n, &k);
    for (int i = 1; i <= k; i++) {
        scanf("%lld", &a[i]); // coeff
        a[i] < 0 && (a[i] += MOD);
    }
    for (int i = 0; i < k; i++) {
        scanf("%lld", &x[i]); // value
        x[i] < 0 && (x[i] += MOD);
    }

    f[k] = 1;
    for (int i = 0; i < k; i++) f[i] = a[k - i] ? MOD - a[k - i] : 0;

    t[1] = 1;
    polyPow(k + 1, t, f, n, m);

    long long ans = 0;
    for (int i = 0; i < k; i++) ans = (ans + m[i] * x[i]) % MOD;
    printf("%lld\n", ans);

    return 0;
}
```
