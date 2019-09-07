# IO 优化

```c++
struct IO {
    static constexpr int BUF_SIZE = 1048576;
    char ibuf[BUF_SIZE], *ip, *ie, obuf[BUF_SIZE], *op;
    std::streambuf *rfb, *ofb;

    IO() : op(obuf) {
        ip = ie = NULL;
        rfb = std::cin.rdbuf();
        ofb = std::cout.rdbuf();
    }
    ~IO() {
        ofb->sputn(obuf, op - obuf);
    }

    char read() {
        if (ip == ie) ie = (ip = ibuf) + rfb->sgetn(ibuf, BUF_SIZE);
        return ip == ie ? -1 : *ip++;
    }

    void read(int &x) {
        static char c;
        static bool isneg;

        for (isneg = false; !std::isdigit(c) && c != -1; c = read()) isneg |= c == '-';
        for (x = 0; std::isdigit(c) && c != -1; x = x * 10 + (c ^ 48), c = read()) {}

        x = (isneg ? -x : x);
    }

    template <typename T>
    IO &operator>>(T &x) {
        read(x);
        return *this;
    }

    void print(char c) {
        if (op == obuf + BUF_SIZE) {
            ofb->sputn(obuf, BUF_SIZE);
            op = obuf;
        }
        *op++ = c;
    }

    template <typename T>
    void print(T x) {
        static char buf[21];
        static int cnt;
        if (x < 0) {
            print('-');
            x = -x;
        }
        cnt = 0;
        do {
            buf[++cnt] = x % 10 | 48;
        } while (x /= 10);
        while (cnt) print((char) buf[cnt--]);
    }

    template <typename T>
    IO &operator<<(const T &x) {
        print(x);
        return *this;
    }
} in, out;
```
