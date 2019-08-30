[TOC]



## 1.Head

```cpp
#include<bits/stdc++.h>
using namespace std;
//#include<ext/rope>
//using namespace __gnu_cxx
//#include<ext/pb_ds/priority_queue.hpp>
//using namespace __gnu_pbds;
#define lowbit(x) ((x)&(-(x)))
#define pb(x) push_back(x)
#define all(x) (x).begin(),(x).end()
#define clr(a,b) memset(a,b,sizeof(a))
#define caze(T) for(cin>>T;T;T--)
#define Endl ('\n')
#define fi first
#define se second
#define ll long long
#define ull unsigned long long
#define i16 __int128
#define pll pair<ll,ll>
#define pli pair<ll,int>
#define pil pair<int,ll>
#define pii pair<int,int>
#define IOS ios::sync_with_stdio(0),cin.tie(0),cout.tie(0)
```

<div style="page-break-after:always;"></div>
## 2.Tarjan

### 边双连通分量

```cpp
//low[u]==low[v]:u与v在同一联通分量;
const int maxn=5007;
struct EDGE{int v,nxt;}edge[1000010];
int _,head[maxn];
void AE(int u,int v){edge[_]={v,head[u]},head[u]=_++;}
int n,m,ttime,sdex,col;
int dfn[maxn],low[maxn],stk[maxn],belong[maxn];
bool vis[maxn];
void ini()
{
	_=ttime=sdex=col=0;
	clr(dfn,0);
	clr(head,-1);
	clr(vis,0);
	clr(belong,0);
}
void tarjan(int u,int fa)
{
	dfn[u]=low[u]=++ttime;
	vis[u]=1;
	stk[++sdex]=u;
	for (int i = head[u]; ~i; i=edge[i].nxt)
	{
		int v=edge[i].v;
		if(v==fa) continue;
		if(!dfn[v])
		{
			tarjan(v,u);
			low[u]=min(low[u],low[v]);
		}
		else if(vis[v])
			low[u]=min(low[u],dfn[v]);
	}
	if (low[u]==dfn[u])
	{
		col++;
		do
		{
			belong[stk[sdex]]=col;
			vis[stk[sdex--]]=0;
		}while(vis[u]);
	}
	return;
}
```

<div style="page-break-after:always;"></div>
### 点双连通分量

```cpp
const int maxn=10007;
struct EDGE{int v,nxt;}edge[1000010];
int _,head[maxn];
void AE(int u,int v){edge[_]={v,head[u]},head[u]=_++;}
int n,m,ttime,idx,col,dfn[maxn],low[maxn],stk[maxn],belong[maxn];
bool vis[maxn],mark[maxn];
void init()
{
	_=ttime=idx=col=0;
	clr(dfn,0);
	clr(head,-1);
	clr(vis,0);
	clr(belong,0);
}
void tarjan(int u,int fa)
{
	dfn[u]=low[u]=++ttime;
	vis[u]=1;
	stk[++idx]=u;
	for(int i=head[u];~i;i=edge[i].nxt)
	{
		int v=edge[i].v;
		if(v==fa) continue;
		if(!dfn[v])
		{
			tarjan(v,u);
			low[u]=min(low[v],low[u]);
			if(low[v]>=dfn[u])				//u是割点
			{
				clr(mark,0);
				int tp,p=0,num=0;
				col++;
				do
				{
					tp=stk[idx--];
					mark[tp]=1;
					vis[tp]=0;
					belong[tp]=col;
				} while (tp!=v);
				mark[u]=1;
			}
		}
		else if(vis[v])
			low[u]=min(low[u],dfn[v]);
	}
}
```

<div style="page-break-after:always;"></div>
## 3.二分匹配

### KM

```cpp
const int N=333;

int n,nx,ny;
int link[N],lx[N],ly[N],slack[N];
int visx[N],visy[N],w[N][N];
bool dfs(int x)
{
	visx[x]=1;
	for(int y=0;y<ny;++y)
	{
		if(visy[y])
			continue;
		int tp=lx[x]+ly[y]-w[x][y];
		if(tp==0)
		{
			visy[y]=1;
			if(link[y]==-1||dfs(link[y]))
			{
				link[y]=x;
				return 1;
			}
		}
		else if(slack[y]>tp)
			slack[y]=tp;
	}
	return 0;
}
int KM()
{
	clr(link,-1);
	clr(ly,0);
	for(int i=0;i<nx;++i)
	{
		lx[i]=-inf;
		for(int j=0;j<ny;++j)
			if(w[i][j]>lx[i])
				lx[i]=w[i][j];
	}
	for(int x=0;x<nx;++x)
	{
		for(int i=0;i<ny;++i)
			slack[i]=inf;
		while(1)
		{
			clr(visx,0);
			clr(visy,0);
			if(dfs(x))
				break;
			int d=inf;
			for(int i=0;i<ny;++i)
				if(!visy[i]&&d>slack[i])
					d=slack[i];
			for(int i=0;i<nx;++i)
				if(visx[i])
					lx[i]-=d;
			for(int i=0;i<ny;++i)
				if(visy[i])
					ly[i]+=d;
				else
					slack[i]-=d;
		}
	}
	int ret=0;
	for(int i=0;i<ny;++i)
		if(link[i]!=-1)
			ret+=w[link[i]][i];
	return ret;
}
void solve(int n)
{
	for(int i=0;i<n;++i)
		for(int j=0;j<n;++j)
			scanf("%d",&w[i][j]);
	nx=ny=n;
	printf("%d\n",KM());
}
```

<div style="page-break-after:always;"></div>
### 匈牙利

```cpp
const int N=255;
vector<int>v[N];
bool used[N];
int lef[N],n;
void init()
{
    for(int i=1;i<=H;++i)
        v[i].clear();
    clr(lef,-1);
}
bool dfs(int x)
{
    for(auto c:v[x])
    {
        if(used[c]==0)
        {
            used[c]=1;
            if(lef[c]==-1||dfs(lef[c]))
            {
                lef[c]=x;
                return 1;
            }
        }
    }
    return 0;
}
int maxmatch()
{
	int ret=0;
	for(int i=1;i<=n;++i)
	{
		clr(used,0);
		ret+=dfs(i);
	}
	return ret;
}
```

<div style="page-break-after:always;"></div>
### HK

```cpp
const int N=3010;
int mx[N],my[N],nx,ny;
int dx[N],dy[N],dis;
bool vis[N];
vector<int>v[N];
int n,m;
bool bfs()
{
	queue<int>Q;
	dis=inf;
	clr(dx,-1);
	clr(dy,-1);
	for(int i=0;i<nx;++i)
		if(mx[i]==-1)
			Q.push(i),dx[i]=0;
	while(!Q.empty())
	{
		int u=Q.front();
		Q.pop();
		if(dx[u]>dis)
			break;
		for(auto w:v[u])
			if(dy[w]==-1)
			{
				dy[w]=dx[u]+1;
				if(my[w]==-1)
					dis=dy[w];
				else
				{
					dx[my[w]]=dy[w]+1;
					Q.push(my[w]);
				}
			}
	}
	return dis!=inf;
}
bool dfs(int u)
{
	for(auto w:v[u])
		if(!vis[w]&&dy[w]==dx[u]+1)
		{
			vis[w]=1;
			if(my[w]!=-1&&dy[w]==dis)
				continue;
			if(my[w]==-1||dfs(my[w]))
			{
				my[w]=u;
				mx[u]=w;
				return 1;
			}
		}
	return 0;
}
int maxmatch()
{
	nx=n,ny=m;
	int ret=0;
	clr(mx,-1);
	clr(my,-1);
	while(bfs())
	{
		clr(vis,0);
		for(int i=0;i<nx;++i)
			if(mx[i]==-1)
				ret+=dfs(i);
	}
	return ret;
}
```

<div style="page-break-after:always;"></div>
## 4.网络流

### 费用流(dij)

```cpp
const int MAXN=111;
struct EDGE{int to,cap,cost,flow,nxt;}edge[1001000];
int head[MAXN];
int _;
void AE(int from,int to,int cap,int cost)
{
	edge[_]={to,cap,cost,0,head[from]},head[from]=_++;
	edge[_]={from,0,-cost,0,head[to]},head[to]=_++;
}
int cost,flow,h[MAXN],dist[MAXN],pre[MAXN],N;
void init(int x)
{
	N=x;
	_=0;
	clr(head,-1);
}
void MCMF(int s,int t,int f)
{
	fill(h,h+N+1,0);
	while(f)
	{
		priority_queue<pii,vector<pii>,greater<pii> >q;
		clr(pre,-1);
		clr(dist,INF);
		dist[s]=0;
		q.push(pii(0,s));
		while(!q.empty())
		{
			pii now=q.top();q.pop();
			if(dist[now.second]<now.first) continue;
			int u=now.second;
			for(int i=head[u];~i;i=edge[i].nxt)
			{
				EDGE &e=edge[i];
				if(e.cap>e.flow&&dist[e.to]>dist[u]+e.cost+h[u]-h[e.to])
				{
					dist[e.to]=dist[u]+e.cost+h[u]-h[e.to];
					pre[e.to]=i;
					q.push({dist[e.to],e.to});
				}
			}
		}
		if(dist[t]==inf) break;
		for (int i = 0; i <= N; ++i)
			h[i]+=dist[i];
		int d=f;
		for (int i = pre[t]; ~i; i=pre[edge[i^1].to])
			d=min(d,edge[i].cap-edge[i].flow);
		f-=d,flow+=d;
		cost+=d*h[t];
		for (int i = pre[t]; ~i; i=pre[edge[i^1].to])
		{
			edge[i].flow+=d;
			edge[i^1].flow-=d;
		}
	}
}
```

<div style="page-break-after:always;"></div>
### 费用流(SPFA)

```cpp
const int MAXN=800010;
struct EDGE{int to,cap,cost,flow,nxt;}edge[MAXN<<2];
int head[MAXN];
int tot;
void AE(int from,int to,int cap,int cost)
{
	edge[tot]={to,cap,cost,0,head[from]},head[from]=tot++;
	edge[tot]={from,0,-cost,0,head[to]},head[to]=tot++;
}
int cost,flow;
int h[MAXN];
int dis[MAXN],pre[MAXN];
int N;
void init(int n)
{
	N=n;
	tot=0;
	memset(head,-1,sizeof(head));
}

/*
 * SPFA 算法判断是否存在s到t的通路
 */
bool spfa(int s,int t)
{
	queue<int>q;
	bool vis[MAXN];
	for(int i=0; i<N; i++)
	{
		dis[i]=inf;
		vis[i]=false;
		pre[i]=-1;
	}
	dis[s]=0;
	vis[s]=true;
	q.push(s);
	while(!q.empty())
	{
		int u=q.front();
		q.pop();
		vis[u]=false;
		for(int i=head[u]; i!=-1; i=edge[i].nxt)   //遍历所有与u临接的点
		{
			int v=edge[i].to;
			if(edge[i].cap>edge[i].flow&&dis[v]>dis[u]+edge[i].cost)    //如果可以松弛该点
			{
				dis[v]=dis[u]+edge[i].cost;
				pre[v]=i;
				if(!vis[v]) //如果该点不在队列中，入队
				{
					vis[v]=true;
					q.push(v);
				}
			}
		}
	}
	return (pre[t]!=-1);    //返回是否s到t是否有路径
}

/*
 * int s 起点
 * int t 终点
 * 返回费用cost
 */
int minCostMaxFlow(int s,int t)
{
	int cost=0;
	while(spfa(s,t))
	{
		int minn=inf;   //当前路径可增广值
		for(int i=pre[t]; i!=-1; i=pre[edge[i^1].to])   //因为建图时每增加一条边会同时加入它的反向边，因此i^1为找出与i刚好相反的部分
		{
			if(minn>edge[i].cap-edge[i].flow)
				minn=edge[i].cap-edge[i].flow;
		}
		for(int i=pre[t]; i!=-1; i=pre[edge[i^1].to])   //修改图，计算花费
		{
			edge[i].flow+=minn;
			edge[i^1].flow-=minn;
			cost+=edge[i].cost*minn;
		}
	}
	return cost;
}
```

<div style="page-break-after:always;"></div>
### 最大流

```cpp
const int N=505;
struct EDGE{int v,cap,nxt;}edge[1000010];
int head[N];int _;
int level[N];
void init(){clr(head,-1);_=0;}
void AE(int u,int v,int w)
{
	edge[_]={v,w,head[u]};head[u]=_++;
	edge[_]={u,0,head[v]};head[v]=_++;
}
bool bfs(int s,int t)
{
	clr(level,-1);
	level[s]=0;
	queue<int>q;
	q.push(s);
	while(!q.empty())
	{
		int u=q.front();q.pop();
		for(int i=head[u];~i;i=edge[i].nxt)
		{
			int v=edge[i].v,f=edge[i].cap;
			if(level[v]==-1&&f)
			{
				level[v]=level[u]+1;
				if(v==t) return 1;
				q.push(v);
			}
		}
	}
	return ~level[t];
}
int dfs(int u,int t,int flow)
{
	if(u==t||!flow) return flow;
	int ret=0;
	for(int i=head[u];~i;i=edge[i].nxt)
	{
		EDGE &e=edge[i];
		int v=e.v;int f=e.cap;
		if((!f)||level[v]!=level[u]+1) continue;
		int fee=dfs(v,t,min(f,flow));
		if(fee<=0) continue;
		edge[i].cap-=fee,edge[i^1].cap+=fee;
		flow-=fee;ret+=fee;
		if(!flow) break;
	}
	if(!ret) level[u]=inf;
	return ret;
}
int Dinic(int s,int t)
{
	int ans=0;
	while(bfs(s,t))
		ans+=dfs(s,t,INT_MAX);
	return ans;
}
```

<div style="page-break-after:always;"></div>
## 5.几何

```cpp
const double pi=acos(-1.0);
const double eps=1e-8;
int dcmp(double x){return fabs(x)<=eps?0:(x<0?-1:1);}
double sqr(double x){return x*x;}
struct point
{
	double x,y,id;
	point(){}
	point(double x,double y,int id=-1):x(x),y(y),id(id) {}
	point operator-(const point w)const {return point(x-w.x,y-w.y);}
	point operator+(const point w)const {return point(x+w.x,y+w.y);}
	double operator*(const point& w)const {return x*w.x+y*w.y;}
	point operator*(double a) {return point(x*a,y*a);}
	double operator^(const point& w)const {return x*w.y-y*w.x;}
	point operator/(double a) {return point(x/a,y/a);}
	friend ostream &operator<<(ostream& out,const point& w) {out<<'('<<w.x<<','<<w.y<<')';return out;}
	void input(){scanf("%lf%lf",&x,&y);}
	double len2(){return x*x+y*y;}
	double len(){return sqrt(x*x+y*y);}
	point change_len(double r)
	{
		double l=len();
		if(dcmp(l)==0) return *this;
		r/=l;
		return point(x*r,y*r);
	}
};
inline double cross(const point& A,const point& B){return A.x*B.y-B.x*A.y;}
inline double dot(const point& q,const point& w){return q.x*w.x+q.y*w.y;}
inline double Xmul(const point& A,const point& B,const point& C){return cross(C-A,B-A);}
inline double dis(const point& q,const point& w){return sqrt(dot(q-w,q-w));}
inline double rad(const point& A,const point& B){return fabs(atan2(fabs(cross(A,B)),dot(A,B)));}
int Andrew(int n,point *st,point *ed)
{
	sort(st,st+n,[](const point& A,const point& B)->bool{return A.x==B.x?A.y<B.y:A.x<B.x;});
	int tot=0;
	for (int i = 0; i < n; ++i)
	{
		while(tot>1&&cross(ed[tot-1]-ed[tot-2],st[i]-ed[tot-2])<=0) tot--;
		ed[tot++]=st[i];
	}
	int tp=tot;
	for (int i = n - 2; ~i; --i)
	{
		while(tot>tp&&cross(ed[tot-1]-ed[tot-2],st[i]-ed[tot-2])<=0) tot--;
		ed[tot++]=st[i];
	}
	tot-=(n>1);
	return tot;
}
double Area(int n,point *p)
{
	double S=0;
	for (int i = 1; i < n - 1; ++i)
		S+=fabs(Xmul(p[0],p[i],p[i+1]));
	return S/2;
}
struct Line
{
	point u,v;
	double k;
	Line(){}
	Line(point u,point v):u(u),v(v){k=atan2(v.y-u.y,v.x-u.x);}
	Line(point u,double k):u(u),k(k){v=u+(dcmp(k-pi/2)?point(1,tan(k)):point(0,1));}
	void input(){u.input();v.input();get_angle();}
	void get_angle(){k=atan2(v.y-u.y,v.x-u.x);}
	double len(){return dis(u,v);}
	double pdis(point w) {return fabs(cross(w-u,v-u)/len());}
	point operator&(const Line& b)const
	{
		point ret=u;
		double t=(cross(u-b.u,b.u-b.v))/cross(u-v,b.u-b.v);
		ret.x+=(v.x-u.x)*t;
		ret.y+=(v.y-u.y)*t;
		return ret;
	}
	point project(const point& w)const{return u+(((v-u)*((v-u)*(w-u)))/(v-u).len2());}
	friend ostream &operator<<(ostream &out,const Line& w){out<<w.u<<"->"<<w.v;return out;}
};
Line Q[100010];
void Hpi(int n,Line *line,point *res,int &resn)
{
	for (int i = 0; i < n; ++i) line[i].get_angle();
	int tot=n;
	sort(line,line+n,[](const Line& A,const Line& B)->bool{return fabs(A.k-B.k)>eps?A.k<B.k:cross(A.u-B.u,B.v-B.u)<0;});
	tot=1;
	for (int i = 1; i < n; ++i)
		if(fabs(line[i].k-line[i-1].k)>eps)
			line[tot++]=line[i];
	int head=0,tail=1;
	Q[0]=line[0];
	Q[1]=line[1];
	resn=0;
	for (int i = 2; i < tot; ++i)
	{
		if(fabs(cross(Q[tail].v-Q[tail].u,Q[tail-1].v-Q[tail-1].u))<eps||fabs(cross(Q[head].v-Q[head].u,Q[head+1].v-Q[head+1].u))<eps)
			return;
		while(head<tail&&(cross((Q[tail]&Q[tail-1])-line[i].u,line[i].v-line[i].u))>eps) tail--;
		while(head<tail&&(cross((Q[head]&Q[head+1])-line[i].u,line[i].v-line[i].u))>eps) head++;
		Q[++tail]=line[i];
	}
	while(head<tail&&(cross(((Q[tail]&Q[tail-1])-Q[head].u),Q[head].v-Q[head].u))>eps) tail--;
	while(head<tail&&(cross(((Q[head]&Q[head-1])-Q[tail].u),Q[tail].v-Q[tail].v))>eps) head++;
	if(tail<=head+1)
		return;
	for (int i = head; i < tail; ++i)
		res[resn++]=Q[i]&Q[i+1];
	if(head<tail-1)
		res[resn++]=Q[head]&Q[tail];
}
struct Circle
{
	point o;
	double r;
	Circle(){}
	Circle(point o,double r):o(o),r(r){}
};
int relation(point w,Line l)
{
	//1:左侧 2:右侧 3:线上
	int c=dcmp(cross(w-l.u,l.v-l.u));
	return c<0?1:(c==0?3:2);
}
int relation(point p,Circle a)
{
	//0:圆外,1:圆上,2:圆内
	double d=dis(p,a.o)-a.r;
	if(dcmp(d)==0) return 1;
	return (dcmp(d)<0?2:0);
}
int relation(Line a,Circle b)
{
	//0:相离,1:相切,2:相交
	double p=a.pdis(b.o);
	if (dcmp (p-b.r) == 0) return 1;
	return (dcmp (p-b.r) < 0 ? 2 : 0);
}
int line_cirlce_intersection(Line l,Circle c,point& p1,point& p2)
{
	if(!relation(l,c))
		return 0;
	point a=l.project(c.o);
	double d=l.pdis(c.o);
	d=sqrt(c.r*c.r-d*d);
	if(dcmp(d)==0)
	{
		p1=a,p2=a;
		return 0;
	}
	p1=a+(l.v-l.u).change_len(d);
	p2=a-(l.v-l.u).change_len(d);
	return 2;
}
double circle_traingle_area(point a,point b,Circle c)
{
	point p=c.o;double r=c.r;
	if(dcmp(cross(p-a,p-b))==0)
		return 0;
	point q[6];
	int len=0;
	q[len++]=a;
	Line l=Line(a,b);
	if (line_cirlce_intersection (l, c, q[1], q[2]) == 2) 
	{
		if (dcmp (dot (a-q[1], b-q[1])) < 0) q[len++] = q[1];
		if (dcmp (dot (a-q[2], b-q[2])) < 0) q[len++] = q[2];
    }
    q[len++]=b;
    if(len==4&&dcmp(dot (q[0]-q[1], q[2]-q[1])) > 0)
    	swap(q[1],q[2]);
    double ans=0;
    for (int i = 0; i < len - 1; ++i)
    {
    	if(relation(q[i],c)==0||relation(q[i+1],c)==0)
    	{
    		double arg=rad(q[i]-p,q[i+1]-p);
    		ans+=r*r*arg/2.0;
    	}
    	else
    		ans+=fabs(cross (q[i]-p, q[i+1]-p))/2;
    }
    return ans;
}
double area_polygon_circle(Circle c,point* p,int n)
{
	double ans=0;
	p[n]=p[0];
	for (int i = 0; i < n; ++i)
	{
		if(dcmp(cross(p[i+1]-c.o,p[i]-c.o))>=0)
			ans+=circle_traingle_area(p[i],p[i+1],c);
		else
			ans-=circle_traingle_area(p[i],p[i+1],c);
	}
	return fabs(ans);
}
```

<div style="page-break-after:always;"></div>
## 6.暴搜

### 第k短路径,不限起点终点,边可经过无数次

```cpp
const int N=50050;
typedef array<ll,4> al4;
vector<pii>g[N];
int n,m,mk;
priority_queue<ll>Q;
priority_queue<al4,vector<al4>,greater<al4> >q;
bool ins(al4 a)
{
	if(Q.size()==mk)
	{
		if(a[0]>=Q.top()) return 0;
		Q.pop();
	}
	Q.push(a[0]);
	return 1;
}
int qq[N];
ll as[N];
int main()
{
	int QAQ,que;caze(QAQ)
	{
		scanf("%d%d%d",&n,&m,&que);
		while(!Q.empty()) Q.pop();
		while(!q.empty()) q.pop();
		for(int i=1;i<=n;++i)
			g[i].clear();
		for(int i=0,u,v,w;i<m;++i)
		{
			scanf("%d%d%d",&u,&v,&w);
			g[u].pb(pii(w,v));
		}
		mk=0;
		for(int i=0;i<que;++i)
		{
			scanf("%d",qq+i);
			mk=max(mk,qq[i]);
		}
		for(int i=1;i<=n;++i)
		{
			sort(all(g[i]));
			if(!g[i].empty())
				q.push({g[i][0].fi,i,g[i][0].se,1});
		}
		while(!q.empty())
		{
			al4 now=q.top();q.pop();
			if(!ins(now)) break;
			int u=now[1],v=now[2],sz=now[3];
			ll val=now[0];
			if(sz<g[u].size())
				q.push({val-g[u][sz-1].fi+g[u][sz].fi,u,g[u][sz].se,sz+1});
			if(!g[v].empty())
				q.push({val+g[v][0].fi,v,g[v][0].se,1});
		}
		for(int i=mk;i;--i)
			as[i]=Q.top(),Q.pop();
		for(int i=0;i<que;++i)
			printf("%lld\n",as[qq[i]]);
	}
}
```

<div style="page-break-after:always;"></div>
### 第k小团 $1\leq N\leq 100,1\leq K\leq 10^6$

```cpp
bitset<111>b[111],now[1000010];
pli p[111];
priority_queue<array<ll,3>,vector<array<ll,3> >,greater<array<ll,3> > >q;
int main()
{
	int n,k;
	scanf("%d%d",&n,&k);
	for(int i=0,x;i<n;++i)
	{
		scanf("%d",&x);
		p[i]=pli(x,i);
	}
	for(int i=0;i<n;++i)
	{
		scanf("%s",s);
		b[i].reset();
		for(int j=0;j<n;++j)
			if(s[j]=='1')
				b[i][j]=1;
	}
	sort(p,p+n);
	if(--k==0) return 0*puts("0");
	int c=0;
	for(int i=0;i<n;++i)
		now[c][i]=1;
	q.push({p[0].fi,c,0});
	while(!q.empty())
	{
		array<ll,3> u=q.top();q.pop();
		ll t=u[0];
		if(--k==0)
			return (printf("%lld\n",t),0);
		int id=u[1],lst=u[2];
		for(int i=lst+1;i<n;++i)
			if(now[id][p[i].se])
			{
				q.push({t-p[lst].fi+p[i].fi,id,i});
				break;
			}
		now[++c]=now[id]&b[p[lst].se];
		for(int i=lst+1;i<n;++i)
			if(now[c][p[i].se])
			{
				q.push({t+p[i].fi,c,i});
				break;
			}
	}
	puts("-1");
}
```

<div style="page-break-after:always;"></div>
### 二分图最大独立集(输出方案)

```cpp
const int maxn=5005;
int mat[maxn];
vector<int>v[maxn];
int vis[maxn];
int n,T;
bool dfs(int x)
{
	vis[x]=T;
	for(int w:v[x])
	{
		if(vis[w]==T) continue;
		vis[w]=T;
		if(mat[w]==-1||dfs(mat[w]))
		{
			mat[w]=x;
			mat[x]=w;
			return 1;
		}
	}
	return 0;
}
void init()
{
	for(int i=0;i<n;++i)
		v[i].clear();
	clr(mat,-1);
}
int a[5005];
bool chk(int x){return (x&(x-1))==0;}
int main()
{
	scanf("%d",&n);
	init();
	for(int i=0;i<n;++i)
	{
		scanf("%d",a+i);
		for(int j=0;j<i;++j)
			if(chk(a[i]^a[j]))
			{
				if(__builtin_parity(a[i]))
					v[i].pb(j);
				else
					v[j].pb(i);
			}
	}
	int ans=0;
	T=0;
	for(int i=0;i<n;++i)
	{
		++T;
		ans+=dfs(i);
	}
	ans=n-ans;
	printf("%d\n",ans);
	++T;
	for(int i=0;i<n;++i)
		if(__builtin_parity(a[i])&&mat[i]==-1)
			dfs(i);
	for(int i=0;i<n;++i)
	{
		if(__builtin_parity(a[i])&&vis[i]==T)
			printf("%d ",a[i]);
		else if(!__builtin_parity(a[i])&&vis[i]!=T)
			printf("%d ",a[i]);
	}
	putchar('\n');
}
```

<div style="page-break-after:always;"></div>
## 7.奇怪的图论问题

### 点数为k最长路径(边带权) $2\leq k\leq 6$

```cpp
//给点染色,dp[S][u]表示状态为S,点u结尾的最长路.
const int N=10010,K=6;
const ll inf=~0ULL>>1;
struct EDGE{int v,w,nxt;}edge[N<<1];
int head[N],_;
void init(){clr(head,-1);_=0;}
void AE(int u,int v,int w){edge[_]={v,w,head[u]};head[u]=_++;}
int c[N],n,m,k,cnt;
ll dp[1<<K][N];
default_random_engine e(time(0));
ll go()
{
	uniform_int_distribution<int>x(0,k-1);
	for(int i=0;i<1<<k;++i)
		for(int j=1;j<=n;++j)
			dp[i][j]=-1;
	for(int i=1;i<=n;++i)
	{
		c[i]=x(e);
		dp[1<<c[i]][i]=0;
	}
	for(int S=1;S<1<<k;++S)
		for(int u=1;u<=n;++u)
		{
			if(((S>>c[u])&1)) continue;
			for(int i=head[u];~i;i=edge[i].nxt)
			{
				int v=edge[i].v,w=edge[i].w;
				if(!((S>>c[v])&1)||(dp[S][v]==-1)) continue;
				dp[S|(1<<c[u])][u]=max(dp[S][v]+w,dp[S|(1<<c[u])][u]);
			}
		}
	ll ret=-1;
	for(int i=1;i<=n;++i)
		ret=max(ret,dp[(1<<k)-1][i]);
	return ret;
}
int main()
{
	int QAQ;caze(QAQ)
	{
		init();
		scanf("%d%d%d",&n,&m,&k);
		for(int i=0;i<m;++i)
		{
			int u,v,w;
			scanf("%d%d%d",&u,&v,&w);
			AE(u,v,w);AE(v,u,w);
		}
		int T=200;
		ll ans=-1;
		for(int i=0;i<T;++i)
			ans=max(ans,go());
		if(ans==-1) puts("impossible");
		else printf("%lld\n",ans);
	}
}
```

<div style="page-break-after:always;"></div>
### 拟图(拟阵交)

```cpp
const int N=88,M=1<<18,inf=~0U>>1;
int h[N],nxt[N<<1],v[N<<1],vis[N],num[N],lim[N],use[N],_,ans;
int from[N],cost[N];
int n,m,cnt;
void AE(int x,int y){v[++_]=y,nxt[_]=h[x],h[x]=_;}
void dfs(int u)
{
	if(vis[u]) return;
	vis[u]=1;
	for(int i=h[u];i;i=nxt[i])
		if(use[i>>1])
			dfs(v[i]);
}
bool chk1()
{
	for(int i=0;i<=n;++i) vis[i]=0;
	dfs(0);
	for(int i=0;i<=n;++i) if(!vis[i]) return 0;
	return 1;
}
bool chk2()
{
	for(int i=1;i<=cnt;++i)
		if(num[i]<lim[i]) return 0;
	return 1;
}
namespace Matroid
{
	int h[N],nxt[M],v[M],q[M],_,hd,tl,d[N],pre[N],w[N],in[N];
	void AE(int x,int y){v[++_]=y,nxt[_]=h[x],h[x]=_;}
	bool find()
	{
		int S=m+1,T=S+1;_=0;
		for(int i=0;i<=T;++i) h[i]=w[i]=0,d[i]=inf,in[i]=0;
		for(int i=1;i<=m;++i)
		{
			w[i]=cost[i];
			if(!use[i]) continue;
			w[i]=-w[i];
			use[i]^=1;num[from[i]]--;
			if(chk1()) AE(S,i);
			if(chk2()) AE(i,T);
			use[i]^=1;num[from[i]]++;
		}
		for(int i=1;i<=m;++i)
			for(int j=1;j<=m;++j)
			{
				if((!use[i])||use[j]) continue;
				use[i]^=1;num[from[i]]--;
				use[j]^=1;num[from[j]]++;
				if(chk1()) AE(j,i);
				if(chk2()) AE(i,j);
				use[i]^=1;num[from[i]]++;
				use[j]^=1;num[from[j]]--;
			}
		hd=tl=0;q[tl++]=S;
		in[S]=1;d[S]=0;
		while(hd<tl)
		{
			int u=q[hd++];
			for(int i=h[u];i;i=nxt[i])
			{
				if(d[v[i]]<=d[u]+w[v[i]]) continue;
				d[v[i]]=d[u]+w[v[i]];
				pre[v[i]]=u;
				if(in[v[i]]) continue;
				q[tl++]=v[i];
				in[v[i]]=1;
			}
			in[u]=0;
		}
		if(d[T]==inf) return 0;
		ans+=d[T];
		while(pre[T]!=S)
		{
			T=pre[T];
			num[from[T]]+=use[T]?-1:1;
			use[T]^=1;
		}
		return 1;
	}
}
int main()
{
	int QAQ;caze(QAQ)
	{
		m=ans=0;_=1;
		int up=0,down=0;
		scanf("%d%d",&n,&cnt);
		for(int i=0;i<=n;++i) h[i]=0;
		for(int i=1;i<=cnt;++i)
		{
			scanf("%d%d",num+i,lim+i);
			down+=lim[i];
			for(int j=0,x,y;j<num[i];++j)
			{
				from[++m]=i;use[m]=1;
				scanf("%d%d%d",&x,&y,cost+m);
				ans+=cost[m];
				AE(x-1,y),AE(y,x-1);
			}
		}
		if(!chk1())
		{
			puts("-1");
			continue;
		}
		up=m;
		while(up>down)
		{
			if(!Matroid::find()) break;
			--up;
		}
		ans=up==down?ans:-1;
		printf("%d\n",ans);
	}
}
```

<div style="page-break-after:always;"></div>
## 8.长链剖分

### k级祖先(在线)

```cpp
const int N=1<<19;
const int L=23;
struct EDGE{int v,nxt;}edge[N<<1];
int head[N],_;
void init(){clr(head,-1);_=0;}
void AE(int u,int v){edge[_]={v,head[u]},head[u]=_++;}
int son[N],dep[N],md[N],top[N],lg[N];
int fa[N][L];
vector<int>up[N],down[N];
void dfs1(int u,int f)
{
	dep[u]=md[u]=dep[f]+1;
	fa[u][0]=f;
	for(int i=1;i<L;++i) fa[u][i]=fa[fa[u][i-1]][i-1];
	for(int i=head[u];~i;i=edge[i].nxt)
	{
		int v=edge[i].v;
		if(v==f) continue;
		dfs1(v,u);
		if(md[v]>md[u])
			son[u]=v,md[u]=md[v];
	}
}
void dfs2(int u,int tp)
{
	top[u]=tp;
	if(son[u]) dfs2(son[u],tp);
	for(int i=head[u];~i;i=edge[i].nxt)
	{
		int v=edge[i].v;
		if(v==fa[u][0]||v==son[u]) continue;
		dfs2(v,v);
	}
}
int query(int u,int k)
{
	if(!k) return u;
	if(k>=dep[u]) return 0;
	u=fa[u][lg[k]],k^=1<<lg[k];
	int t=dep[u]-dep[top[u]];
	if(t>k) return down[top[u]][t-k];
	else return up[top[u]][k-t];
}
int main()
{
	init();
	lg[0]=-1;
	for(int i=1;i<N;++i) lg[i]=lg[i>>1]+1;
	int n,m,ans=0,u,v;
	scanf("%d",&n);
	for(int i=1;i<n;++i)
	{
		scanf("%d%d",&u,&v);
		AE(u,v);AE(v,u);
	}
	dfs1(1,0);dfs2(1,1);
	for(int i=1;i<=n;++i)
	{
		if(i!=top[i]) continue;
		int len=md[i]-dep[i];
		for(int j=0,k=i;j<=len;++j)
			up[i].pb(k),k=fa[k][0];
		for(int j=0,k=i;j<=len;++j)
			down[i].pb(k),k=son[k];
	}
	scanf("%d",&m);
	for(int i=0;i<m;++i)
	{
		scanf("%d%d",&u,&v);
		u^=ans,v^=ans;
		printf("%d\n",(ans=query(u,v)));
	}
}
```

<div style="page-break-after:always;"></div>
### 树形dp

#### Hotel

```cpp
const int N=1<<18;
const int L=21;
struct EDGE{int v,nxt;}edge[N<<1];
int head[N],_;
void init(){clr(head,-1);_=0;}
void AE(int u,int v){edge[_]={v,head[u]},head[u]=_++;}
int md[N],dep[N],son[N],p[N],len[N];ll *f[N],*g[N],u[N<<2],*e=u+1;
ll ans;
void dfs1(int u,int fa)
{
	p[u]=fa;
	dep[u]=md[u]=dep[fa]+1;
	for(int i=head[u];~i;i=edge[i].nxt)
	{
		int v=edge[i].v;
		if(v==fa) continue;
		dfs1(v,u);
		if(md[v]>md[u])
			md[u]=md[v],son[u]=v;
	}
	len[u]=md[u]-dep[u]+1;
}
void dfs2(int u)
{
	f[u][0]=1;
	if(son[u])
	{
		f[son[u]]=f[u]+1;
		g[son[u]]=g[u]-1;
		dfs2(son[u]);
	}
	for(int i=head[u];~i;i=edge[i].nxt)
	{
		int v=edge[i].v;
		if(v==p[u]||v==son[u]) continue;
		int m=len[v];
		f[v]=e;e+=m<<1;
		g[v]=e;e+=m<<1;
		dfs2(v);
		g[u][0]+=g[v][1];
		for(int j=1;j<=m;++j)
		{
			ans+=f[v][j-1]*g[u][j]+f[u][j]*g[v][j+1];
			g[u][j]+=f[u][j]*f[v][j-1];
			if(j<=m-2) g[u][j]+=g[v][j+1];
		}
		for(int j=1;j<=m;++j)
			f[u][j]+=f[v][j-1];
	}
	ans+=g[u][0];
}
int main()
{
	init();
	int u,v,n;
	ans=0;
	scanf("%d",&n);
	for(int i=1;i<n;++i)
	{
		scanf("%d%d",&u,&v);
		AE(u,v);AE(v,u);
	}
	dfs1(1,0);
	f[1]=e;e+=len[1]<<1;
	g[1]=e;e+=len[1]<<1;
	dfs2(1);
	printf("%lld\n",ans);
}
```

<div style="page-break-after:always;"></div>
#### codeforces1009F

```cpp
const int N=1<<20;
const int L=21;
struct EDGE{int v,nxt;}edge[N<<1];
int head[N],_;
void init(){clr(head,-1);_=0;}
void AE(int u,int v){edge[_]={v,head[u]},head[u]=_++;}
int dep[N],son[N],md[N],pre[N],u[N],*e=u+1,*f[N];
int ans[N];
void dfs1(int u,int fa)
{
	pre[u]=fa;
	md[u]=dep[u]=dep[fa]+1;
	for(int i=head[u];~i;i=edge[i].nxt)
	{
		int v=edge[i].v;
		if(v==fa) continue;
		dfs1(v,u);
		if(md[v]>md[u])
			son[u]=v,md[u]=md[v];
	}
}
void dfs2(int u)
{
	f[u][0]=1;
	if(son[u])
	{
		f[son[u]]=f[u]+1;
		dfs2(son[u]);
		ans[u]=ans[son[u]]+1;
	}
	for(int i=head[u];~i;i=edge[i].nxt)
	{
		int v=edge[i].v;
		if(v==pre[u]||v==son[u]) continue;
		int m=md[v]-dep[v]+1;
		f[v]=e,e+=m;
		dfs2(v);
		for(int j=1;j<=m;++j)
		{
			f[u][j]+=f[v][j-1];
			if(f[u][j]>f[u][ans[u]]||f[u][j]==f[u][ans[u]]&&j<ans[u])
				ans[u]=j;
		}
	}
	if(f[u][ans[u]]==1) ans[u]=0;
}
int main()
{
	init();
	int u,v,n;
	scanf("%d",&n);
	for(int i=1,u,v;i<n;++i)
	{
		scanf("%d%d",&u,&v);
		AE(u,v);AE(v,u);
	}
	dfs1(1,0);
	f[1]=e;
	e+=md[1];
	dfs2(1);
	for(int i=1;i<=n;++i)
		printf("%d\n",ans[i]);
}
```

