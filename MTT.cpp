// luogu-judger-enable-o2
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef double db;
const int N = 1 << 18 | 7, L = 19;
const int M = (1 << 15) - 1;
const int all = 1 << 18;
const db pi = acos(-1.0);
namespace FASTIO {
    inline char nc() {
        static char buf[N], *p1 = buf, *p2 = buf;
        return p1==p2&&(p2=(p1=buf)+fread(buf,1,N,stdin),p1==p2)?EOF:*p1++;
    }
    template<class T>inline int scan(T &x) {
        char c; int sign = 1;
        while(!isdigit(c=nc())&&~c)
            if(c=='-') sign *= -1;
        if(!~c) return EOF;
        for(x = 0; isdigit(c); c = nc())
            x = (x<<3) + (x<<1) + (c&15);
        return 1;
    }
}
using namespace FASTIO;
struct cpd {
    db x, y;
    cpd(db _x = 0, db _y = 0):x(_x), y(_y){}
    cpd &operator *= (cpd &R) { return *this = cpd(x*R.x-y*R.y, x*R.y+y*R.x);}
    cpd &operator /= (db t) { return *this = cpd(x/t, y/t); }
    cpd operator + (cpd R) { return cpd(x + R.x, y + R.y); }
    cpd operator - (cpd R) { return cpd(x - R.x, y - R.y); }
    cpd operator * (cpd R) { cpd ret = *this; return ret *= R; }
    friend cpd conj(cpd R) { return cpd(R.x, -R.y); }
};

cpd x[N], y[N], w[N];
void fft(cpd *y, int len, int on) {
    static cpd u, v, t;
    static int r[N], nl;
    int i, j, k, l = __builtin_ctz(len);
    if(nl != l) for(nl = l, i = 0; i < len; i++)
        r[i] = (r[i>>1]>>1)|((i&1)<<l-1);
    for(i = 0; i < len; i++) if(i<r[i]) swap(y[i], y[r[i]]);
    for(i = l = 1; i < len; i <<= 1, l++) {
        for(j = 0; j < len; j += i << 1) {
            for(k = 0; k < i; k++) {
                t = on?conj(w[k<<L-l]):w[k<<L-l];
                u = y[k+j], v = y[k+i+j] * t;
                y[k+j] = u + v, y[k+i+j] = u - v;
            }
        }
    }
    if(on) for(i = 0; i < len; i++) y[i] /= len;
}

int f[N], g[N], h[N], n, m, p;
void solve(int *a, int *b, int *c, int l1, int l2) {
    static cpd dft[4][N];
    static int i, j, d[4];
    int len = l1 + l2 - 1;
    while(len & (len - 1)) len += len & -len;
    for(i = 0; i < len; i++) {
        x[i] = cpd(a[i] & M, a[i] >> 15);
        y[i] = cpd(b[i] & M, b[i] >> 15);
    }
    fft(x, len, 0), fft(y, len, 0);
    const cpd half = cpd(0.5, 0);
    for(i = j = 0; i < len; i++, j = len - i) {
        dft[0][i] = x[i] * y[i] * half;
        dft[1][i] = conj(x[j]) * y[i] * half;
        dft[2][i] = x[i] * conj(y[j]) * half;
        dft[3][i] = conj(x[j]) * conj(y[j]) * half;
    }
    for(i = 0; i < len; i++) {
        x[i] = dft[0][i] + dft[1][i];
        y[i] = dft[2][i] - dft[3][i];
    }
    fft(x, len, 1), fft(y, len, 1);

    for(i = 0; i < len; i++) {
        d[0] = ll(x[i].x + 0.5) % p, d[1] = ll(x[i].y + 0.5) % p;
        d[2] = ll(y[i].y + 0.5) % p, d[3] = ll(y[i].x + 0.5) % p;
        c[i] = ((1ll<<30) * d[3] % p + d[0]) % p;
        c[i] = (c[i] + (1ll<<15) * (d[1] + d[2]) % p) % p;
    }
}

int main() {
#ifdef local
    freopen("in.txt", "r", stdin);
#endif
    for(int i = 0; i <= all; i++)
        w[i] = cpd(cos(i*pi/all), sin(i*pi/all));
    scan(n), scan(m), scan(p);
    n++, m++;
    for(int i = 0; i < n; i++) scan(f[i]);
    for(int i = 0; i < m; i++) scan(g[i]);
    solve(f, g, h, n, m);
    for(int i = 0; i < n + m - 1; i++)
        printf("%d%c", h[i], " \n"[i == n + m - 2]);
    return 0;
}
