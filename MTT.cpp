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
    cpd(db _x=0, db _y=0):x(_x),y(_y){}
    cpd &operator *= (cpd &R) { return *this = cpd{x*R.x-y*R.y, x*R.y+y*R.x};}
    cpd &operator /= (db t) { return *this = cpd{x/t, y/t}; }
    cpd operator + (cpd R) { return cpd{x + R.x, y + R.y}; }
    cpd operator - (cpd R) { return cpd{x - R.x, y - R.y}; }
    cpd operator * (cpd R) { cpd ret = *this; return ret *= R; }
    friend cpd conj(cpd R) { return cpd{R.x, -R.y}; }
    bool operator == (cpd R) { return x == R.x && y == R.y;}
};
cpd x[N], y[N], w[N];

void fft(cpd *y, int len, int on) {
    static int r[N], nl;
    static cpd u, v, t;
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
int t[4][N], res[4][N], f[N], g[N], h[N], n, m, p;
void conv(int *a, int *b, int *c, int len) {
    int i, j;
    for(i = 0; i < len; i++) x[i] = cpd{a[i], b[i]};
    fft(x, len, 0);
    for(i = j = 0; i < len; i++, j = len - i)
        y[i] = (x[i] * x[i] - conj(x[j] * x[j])) * cpd{0, -0.25};
    fft(y, len, 1);
    for(i = 0; i < len; i++) c[i] = ll(y[i].x + 0.5) % p;
}

void solve(int *a, int *b, int *c, int l1, int l2) {
    int len = l1 + l2 - 1;
    while(len & (len - 1)) len += len & -len;
    for(int i = 0; i < len; i++) {
        t[0][i] = a[i] & M, t[1][i] = a[i] >> 15;
        t[2][i] = b[i] & M, t[3][i] = b[i] >> 15;
    }
    conv(t[0], t[2], res[0], len), conv(t[0], t[3], res[1], len);
    conv(t[1], t[2], res[2], len), conv(t[1], t[3], res[3], len);
    for(int i = 0; i < len; i++) {
        c[i] = ((1ll<<30) * res[3][i] % p + res[0][i]) % p;
        (c[i] += (1ll<<15) * (res[1][i] + res[2][i]) % p) %= p;
    }
}

int main() {
#ifdef local
    freopen("in.txt", "r", stdin);
#endif
    for(int i = 0; i <= all; i++)
        w[i] = cpd{cos(i*pi/all), sin(i*pi/all)};
    scan(n), scan(m), scan(p);
    n++, m++;
    for(int i = 0; i < n; i++) scan(f[i]);
    for(int i = 0; i < m; i++) scan(g[i]);
    solve(f, g, h, n, m);
    for(int i = 0; i < n + m - 1; i++)
        printf("%d%c", h[i], " \n"[i == n + m - 2]);
    return 0;
}
