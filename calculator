#include<bits/stdc++.h>
using namespace std;
const int maxn = 1e5 + 7;

char suffix[maxn], mid[maxn];
int prio[1 << 8];

inline bool blank(const char &c) {
	return c == '\n' || c == '\r' || c == ' ' || c == '\t';
}

void trans(const char *c) {
	int cnt = 0, len = strlen(c), t = 0;
	prio['-'] = prio['+'] = 1;
	prio['*'] = prio['/'] = 2;
	stack<char> op;
	for(int i = 0; i < len; i++) {
		if(blank(c[i])) continue;
		if(isdigit(c[i])) {
			suffix[t++] = c[i];
			if(!isdigit(c[i + 1])) suffix[t++] = ' ';
		} else {
			if(c[i] == ')') {
				while(op.top() != '(') suffix[t++] = op.top(), op.pop();
				op.pop();
			} else if(c[i] == '(' || op.empty() || prio[c[i]] > prio[op.top()]) {
				op.push(c[i]);
			} else {
				while(!op.empty() && prio[c[i]] <= prio[op.top()])
					suffix[t++] = op.top(), op.pop();
				op.push(c[i]);
			}
		}
	}
	while(!op.empty()) suffix[t++] = op.top(), op.pop();
	suffix[t] = '\0';
}

int calc(const char *c) {
	int len = strlen(c), tmpnum = 0;
	stack<int> num;
	for(int i = 0; i < len; i++) {
		if(blank(c[i])) continue;
		if(isdigit(c[i])) {
			(tmpnum *= 10) += c[i] - '0';
			if(!isdigit(c[i + 1])) num.push(tmpnum), tmpnum = 0;
		} else {
			int t2 = num.top(); num.pop();
			int t1 = num.top(); num.pop();
			if(c[i] == '*') t1 *= t2;
			else if(c[i] == '/') t1 /= t2;
			else if(c[i] == '+') t1 += t2;
			else if(c[i] == '-') t1 -= t2;
			num.push(t1);
		}
	}
	return num.top();
}

int main() {
	while(~scanf("%s", mid)) {
		trans(mid);
		printf("%d\n", calc(suffix));
	}
	return 0;
}
