#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <vector>
#include <climits>
using namespace std;
 
using ll = long long;
using db = long double; // or double, if TL is tight
using str = string; // yay python! 

// pairs
using pi = pair<int,int>;
using pl = pair<ll,ll>;
using pd = pair<db,db>;
#define mp make_pair
#define f first
#define s second

#define tcT template<class T
#define tcTU tcT, class U
// ^ lol this makes everything look weird but I'll try it
tcT> using V = vector<T>; 
tcT, size_t SZ> using AR = array<T,SZ>; 
using vi = V<int>;
using vb = V<bool>;
using vl = V<ll>;
using vd = V<db>;
using vs = V<str>;
using vpi = V<pi>;
using vpl = V<pl>;
using vpd = V<pd>;

// vectors
// oops size(x), rbegin(x), rend(x) need C++17
#define sz(x) int((x).size())
#define bg(x) begin(x)
#define all(x) bg(x), end(x)
#define rall(x) x.rbegin(), x.rend() 
#define sor(x) sort(all(x)) 
#define rsz resize
#define ins insert 
#define pb push_back
#define eb emplace_back
#define ft front()
#define bk back()

#define lb lower_bound
#define ub upper_bound
tcT> int lwb(V<T>& a, const T& b) { return int(lb(all(a),b)-bg(a)); }
tcT> int upb(V<T>& a, const T& b) { return int(ub(all(a),b)-bg(a)); }

// loops
#define FOR(i,a,b) for (int i = (a); i < (b); ++i)
#define F0R(i,a) FOR(i,0,a)
#define ROF(i,a,b) for (int i = (b)-1; i >= (a); --i)
#define R0F(i,a) ROF(i,0,a)
#define rep(a) F0R(_,a)
#define each(a,x) for (auto& a: x)

const int MOD = (int)1e9+7; // 998244353;
const int MX = (int)2e5+5;
const ll BIG = 1e18; // not too close to LLONG_MAX
const db PI = acos((db)-1);
const int dx[4]{1,0,-1,0}, dy[4]{0,1,0,-1}; // for every grid problem!!
mt19937 rng((uint32_t)chrono::steady_clock::now().time_since_epoch().count()); 
template<class T> using pqg = priority_queue<T,vector<T>,greater<T>>;

// bitwise ops
// also see https://gcc.gnu.org/onlinedocs/gcc/Other-Builtins.html
constexpr int pct(int x) { return __builtin_popcount(x); } // # of bits set
constexpr int bits(int x) { // assert(x >= 0); // make C++11 compatible until USACO updates ...
	return x == 0 ? 0 : 31-__builtin_clz(x); } // floor(log2(x)) 
constexpr int p2(int x) { return 1<<x; }
constexpr int msk2(int x) { return p2(x)-1; }

ll cdiv(ll a, ll b) { return a/b+((a^b)>0&&a%b); } // divide a by b rounded up
ll fdiv(ll a, ll b) { return a/b-((a^b)<0&&a%b); } // divide a by b rounded down

tcT> bool ckmin(T& a, const T& b) {
	return b < a ? a = b, 1 : 0; } // set a = min(a,b)
tcT> bool ckmax(T& a, const T& b) {
	return a < b ? a = b, 1 : 0; } // set a = max(a,b)

tcTU> T fstTrue(T lo, T hi, U f) {
	++hi; assert(lo <= hi); // assuming f is increasing
	while (lo < hi) { // find first index such that f is true 
		T mid = lo+(hi-lo)/2;
		f(mid) ? hi = mid : lo = mid+1; 
	} 
	return lo;
}
tcTU> T lstTrue(T lo, T hi, U f) {
	--lo; assert(lo <= hi); // assuming f is decreasing
	while (lo < hi) { // find first index such that f is true 
		T mid = lo+(hi-lo+1)/2;
		f(mid) ? lo = mid : hi = mid-1;
	} 
	return lo;
}
tcT> void remDup(vector<T>& v) { // sort and remove duplicates
	sort(all(v)); v.erase(unique(all(v)),end(v)); }
tcTU> void erase(T& t, const U& u) { // don't erase
	auto it = t.find(u); assert(it != end(t));
	t.erase(it); } // element that doesn't exist from (multi)set

#define tcTUU tcT, class ...U

inline namespace Helpers {
	//////////// is_iterable
	// https://stackoverflow.com/questions/13830158/check-if-a-variable-type-is-iterable
	// this gets used only when we can call begin() and end() on that type
	tcT, class = void> struct is_iterable : false_type {};
	tcT> struct is_iterable<T, void_t<decltype(begin(declval<T>())),
	                                  decltype(end(declval<T>()))
	                                 >
	                       > : true_type {};
	tcT> constexpr bool is_iterable_v = is_iterable<T>::value;

	//////////// is_readable
	tcT, class = void> struct is_readable : false_type {};
	tcT> struct is_readable<T,
	        typename std::enable_if_t<
	            is_same_v<decltype(cin >> declval<T&>()), istream&>
	        >
	    > : true_type {};
	tcT> constexpr bool is_readable_v = is_readable<T>::value;

	//////////// is_printable
	// // https://nafe.es/posts/2020-02-29-is-printable/
	tcT, class = void> struct is_printable : false_type {};
	tcT> struct is_printable<T,
	        typename std::enable_if_t<
	            is_same_v<decltype(cout << declval<T>()), ostream&>
	        >
	    > : true_type {};
	tcT> constexpr bool is_printable_v = is_printable<T>::value;
}

inline namespace Input {
	tcT> constexpr bool needs_input_v = !is_readable_v<T> && is_iterable_v<T>;
	tcTUU> void re(T& t, U&... u);
	tcTU> void re(pair<T,U>& p); // pairs

	// re: read
	tcT> typename enable_if<is_readable_v<T>,void>::type re(T& x) { cin >> x; } // default
	tcT> void re(complex<T>& c) { T a,b; re(a,b); c = {a,b}; } // complex
	tcT> typename enable_if<needs_input_v<T>,void>::type re(T& i); // ex. vectors, arrays
	tcTU> void re(pair<T,U>& p) { re(p.f,p.s); }
	tcT> typename enable_if<needs_input_v<T>,void>::type re(T& i) {
		each(x,i) re(x); }
	tcTUU> void re(T& t, U&... u) { re(t); re(u...); } // read multiple

	// rv: resize and read vectors
	void rv(size_t) {}
	tcTUU> void rv(size_t N, V<T>& t, U&... u);
	template<class...U> void rv(size_t, size_t N2, U&... u);
	tcTUU> void rv(size_t N, V<T>& t, U&... u) {
		t.rsz(N); re(t);
		rv(N,u...); }
	template<class...U> void rv(size_t, size_t N2, U&... u) {
		rv(N2,u...); }

	// dumb shortcuts to read in ints
	void decrement() {} // subtract one from each
	tcTUU> void decrement(T& t, U&... u) { --t; decrement(u...); }
	#define ints(...) int __VA_ARGS__; re(__VA_ARGS__);
	#define int1(...) ints(__VA_ARGS__); decrement(__VA_ARGS__);
}

inline namespace ToString {
	tcT> constexpr bool needs_output_v = !is_printable_v<T> && is_iterable_v<T>;

	// ts: string representation to print
	tcT> typename enable_if<is_printable_v<T>,str>::type ts(T v) {
		stringstream ss; ss << fixed << setprecision(15) << v;
		return ss.str(); } // default
	tcT> str bit_vec(T t) { // bit vector to string
		str res = "{"; F0R(i,sz(t)) res += ts(t[i]);
		res += "}"; return res; }
	str ts(V<bool> v) { return bit_vec(v); }
	template<size_t SZ> str ts(bitset<SZ> b) { return bit_vec(b); } // bit vector
	tcTU> str ts(pair<T,U> p); // pairs
	tcT> typename enable_if<needs_output_v<T>,str>::type ts(T v); // vectors, arrays
	tcTU> str ts(pair<T,U> p) { return "("+ts(p.f)+", "+ts(p.s)+")"; }
	tcT> typename enable_if<is_iterable_v<T>,str>::type ts_sep(T v, str sep) {
		// convert container to string w/ separator sep
		bool fst = 1; str res = "";
		for (const auto& x: v) {
			if (!fst) res += sep;
			fst = 0; res += ts(x);
		}
		return res;
	}
	tcT> typename enable_if<needs_output_v<T>,str>::type ts(T v) {
		return "{"+ts_sep(v,", ")+"}"; }

	// for nested DS
	template<int, class T> typename enable_if<!needs_output_v<T>,vs>::type 
	  ts_lev(const T& v) { return {ts(v)}; }
	template<int lev, class T> typename enable_if<needs_output_v<T>,vs>::type 
	  ts_lev(const T& v) {
		if (lev == 0 || !sz(v)) return {ts(v)};
		vs res;
		for (const auto& t: v) {
			if (sz(res)) res.bk += ",";
			vs tmp = ts_lev<lev-1>(t);
			res.ins(end(res),all(tmp));
		}
		F0R(i,sz(res)) {
			str bef = " "; if (i == 0) bef = "{";
			res[i] = bef+res[i];
		}
		res.bk += "}";
		return res;
	}
}

inline namespace Output {
	template<class T> void pr_sep(ostream& os, str, const T& t) { os << ts(t); }
	template<class T, class... U> void pr_sep(ostream& os, str sep, const T& t, const U&... u) {
		pr_sep(os,sep,t); os << sep; pr_sep(os,sep,u...); }
	// print w/ no spaces
	template<class ...T> void pr(const T&... t) { pr_sep(cout,"",t...); } 
	// print w/ spaces, end with newline
	void ps() { cout << "\n"; }
	template<class ...T> void ps(const T&... t) { pr_sep(cout," ",t...); ps(); } 
	// debug to cerr
	template<class ...T> void dbg_out(const T&... t) {
		pr_sep(cerr," | ",t...); cerr << endl; }
	void loc_info(int line, str names) {
		cerr << "Line(" << line << ") -> [" << names << "]: "; }
	template<int lev, class T> void dbgl_out(const T& t) {
		cerr << "\n\n" << ts_sep(ts_lev<lev>(t),"\n") << "\n" << endl; }
	#ifdef LOCAL
		#define dbg(...) loc_info(__LINE__,#__VA_ARGS__), dbg_out(__VA_ARGS__)
		#define dbgl(lev,x) loc_info(__LINE__,#x), dbgl_out<lev>(x)
	#else // don't actually submit with this
		#define dbg(...)
		#define dbgl(lev,x)
	#endif

	// https://stackoverflow.com/questions/47980498/accurate-c-c-clock-on-a-multi-core-processor-with-auto-overclock?noredirect=1&lq=1
	const auto beg = std::chrono::high_resolution_clock::now();
	void dbg_time() {
		auto duration = chrono::duration<double>(
			std::chrono::high_resolution_clock::now() - beg);
		dbg(duration.count());
	}
}

inline namespace FileIO {
	void setIn(str s)  { freopen(s.c_str(),"r",stdin); }
	void setOut(str s) { freopen(s.c_str(),"w",stdout); }
	void setIO(str s = "") {
		cin.tie(0)->sync_with_stdio(0); // unsync C / C++ I/O streams
		// cin.exceptions(cin.failbit);
		// throws exception when do smth illegal
		// ex. try to read letter into int
		if (sz(s)) setIn(s+".in"), setOut(s+".out"); // for old USACO
	}
}

/**
 * Description: Use in place of \texttt{complex<T>}.
 * Source: http://codeforces.com/blog/entry/22175, KACTL
 * Verification: various
 */

using T = db; // or ll
const T EPS = 1e-9; // adjust as needed
using P = pair<T,T>; using vP = V<P>; using Line = pair<P,P>;
int sgn(T a) { return (a>EPS)-(a<-EPS); }
T sq(T a) { return a*a; }

bool close(const P& a, const P& b) { 
	return sgn(a.f-b.f) == 0 && sgn(a.s-b.s) == 0; } 
T norm(const P& p) { return sq(p.f)+sq(p.s); }
T abs(const P& p) { return sqrt(norm(p)); }
T arg(const P& p) { return atan2(p.s,p.f); }
P conj(const P& p) { return P(p.f,-p.s); }
P perp(const P& p) { return P(-p.s,p.f); }
P dir(T ang) { return P(cos(ang),sin(ang)); }

P operator-(const P& l) { return P(-l.f,-l.s); }
P operator+(const P& l, const P& r) { 
	return P(l.f+r.f,l.s+r.s); }
P operator-(const P& l, const P& r) { 
	return P(l.f-r.f,l.s-r.s); }
P operator*(const P& l, const T& r) { 
	return P(l.f*r,l.s*r); }
P operator*(const T& l, const P& r) { return r*l; }
P operator/(const P& l, const T& r) { 
	return P(l.f/r,l.s/r); }
P operator*(const P& l, const P& r) { 
	return P(l.f*r.f-l.s*r.s,l.s*r.f+l.f*r.s); }
P operator/(const P& l, const P& r) { 
	return l*conj(r)/norm(r); }
P& operator+=(P& l, const P& r) { return l = l+r; }
P& operator-=(P& l, const P& r) { return l = l-r; }
P& operator*=(P& l, const T& r) { return l = l*r; }
P& operator/=(P& l, const T& r) { return l = l/r; }
P& operator*=(P& l, const P& r) { return l = l*r; }
P& operator/=(P& l, const P& r) { return l = l/r; }

P unit(const P& p) { return p/abs(p); }
T dot(const P& a, const P& b) { return a.f*b.f+a.s*b.s; }
T dot(const P& p, const P& a, const P& b) { return dot(a-p,b-p); }
T cross(const P& a, const P& b) { return a.f*b.s-a.s*b.f; }
T cross(const P& p, const P& a, const P& b) {
	return cross(a-p,b-p); }
P reflect(const P& p, const Line& l) {
	P a = l.f, d = l.s-l.f;
	return a+conj((p-a)/d)*d; }
P foot(const P& p, const Line& l) {
	return (p+reflect(p,l))/(T)2; }
bool onSeg(const P& p, const Line& l) {
	return sgn(cross(l.f,l.s,p)) == 0 && sgn(dot(p,l.f,l.s)) <= 0; }

/**
 * Description: Basic 3D geometry. 
 * Source: Own
 * Verification: (haven't done much 3D geo yet)
	* AMPPZ 2011 Cross Spider
	* https://atcoder.jp/contests/JAG2013Spring/tasks/icpc2013spring_h
	* https://codeforces.com/gym/102040 - I
	* https://codeforces.com/gym/102452/problem/F
 */

// #include "../Primitives/Point.h"

/**
using T = db;
int sgn(T x) { return (x>0)-(x<0); }
T sq(T x) { return x*x; }
*/

using P3 = AR<T,3>; using Tri = AR<P3,3>; using vP3 = V<P3>;
T norm(const P3& x) { 
	T sum = 0; F0R(i,3) sum += sq(x[i]);
	return sum; }
T abs(const P3& x) { return sqrt(norm(x)); }

P3& operator+=(P3& l, const P3& r) { F0R(i,3) l[i] += r[i]; 
	return l; }
P3& operator-=(P3& l, const P3& r) { F0R(i,3) l[i] -= r[i]; 
	return l; }
P3& operator*=(P3& l, const T& r) { F0R(i,3) l[i] *= r; 
	return l; }
P3& operator/=(P3& l, const T& r) { F0R(i,3) l[i] /= r; 
	return l; }
P3 operator-(P3 l) { l *= -1; return l; }
P3 operator+(P3 l, const P3& r) { return l += r; }
P3 operator-(P3 l, const P3& r) { return l -= r; }
P3 operator*(P3 l, const T& r) { return l *= r; }
P3 operator*(const T& r, const P3& l) { return l*r; }
P3 operator/(P3 l, const T& r) { return l /= r; }

P3 unit(const P3& x) { return x/abs(x); }
T dot(const P3& a, const P3& b) { 
	T sum = 0; F0R(i,3) sum += a[i]*b[i]; 
	return sum; }
P3 cross(const P3& a, const P3& b) {
	return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],
			a[0]*b[1]-a[1]*b[0]}; }
P3 cross(const P3& a, const P3& b, const P3& c) {
	return cross(b-a,c-a); }
P3 perp(const P3& a, const P3& b, const P3& c) {
	return unit(cross(a,b,c)); }

bool isMult(const P3& a, const P3& b) { // for long longs
	P3 c = cross(a,b); F0R(i,sz(c)) if (c[i] != 0) return 0; 
	return 1; }
bool collinear(const P3& a, const P3& b, const P3& c) { 
	return isMult(b-a,c-a); }

T DC(const P3&a,const P3&b,const P3&c,const P3&p) { 
	return dot(cross(a,b,c),p-a); }
bool coplanar(const P3&a,const P3&b,const P3&c,const P3&p) { 
	return DC(a,b,c,p) == 0; }
bool op(const P3& a, const P3& b) { 
	int ind = 0; // going in opposite directions?
	FOR(i,1,3) if (std::abs(a[i]*b[i])>std::abs(a[ind]*b[ind])) 
		ind = i;
	return a[ind]*b[ind] < 0;
}
// coplanar points, b0 and b1 on opposite sides of a0-a1?
bool opSide(const P3&a,const P3&b,const P3&c,const P3&d) { 
	return op(cross(a,b,c),cross(a,b,d)); }
// coplanar points, is a in Triangle b
bool inTri(const P3& a, const Tri& b) { 
	F0R(i,3)if(opSide(b[i],b[(i+1)%3],b[(i+2)%3],a))return 0;
	return 1; }

// point-seg dist
T psDist(const P3&p,const P3&a,const P3&b) { 
	if (dot(a-p,a-b) <= 0) return abs(a-p);
	if (dot(b-p,b-a) <= 0) return abs(b-p);
	return abs(cross(p,a,b))/abs(a-b);
}
// projection onto line
P3 foot(const P3& p, const P3& a, const P3& b) { 
	P3 d = unit(b-a); return a+dot(p-a,d)*d; }
// rotate p about axis
P3 rotAxis(const P3& p, const P3& a, const P3& b, T theta) {
	P3 dz = unit(b-a), f = foot(p,a,b); 
	P3 dx = p-f, dy = cross(dz,dx);
	return f+cos(theta)*dx+sin(theta)*dy;
}
// projection onto plane
P3 foot(const P3& a, const Tri& b) {
	P3 c = perp(b[0],b[1],b[2]);
	return a-c*(dot(a,c)-dot(b[0],c)); }
// line-plane intersection
P3 lpIntersect(const P3&a0,const P3&a1,const Tri&b) { 
	P3 c = unit(cross(b[2]-b[0],b[1]-b[0]));
	T x = dot(a0,c)-dot(b[0],c), y = dot(a1,c)-dot(b[0],c);
	return (y*a0-x*a1)/(y-x);
}

/**
 * Description: Incremental 3D convex hull where not all points 
 	* are coplanar. Normals to returned faces point outwards. 
 	* If coordinates are ints at most $B$ then \texttt{T} 
 	* should be large enough to support ints on the order 
 	* of $B^3$. Changes order of points. The number of returned faces
 	* may depend on the random seed, because points that are on the
 	* boundary of the convex hull may or may not be included
 	* in the output.
 * Time: O(N^2), O(N\log N)
 * Source: 
 	* KACTL
 	* https://codeforces.com/blog/entry/73366?#comment-575862 (mango_lassi)
 	* https://codeforces.com/blog/entry/81768 (Monogon)
 	* https://people.csail.mit.edu/indyk/6.838-old/handouts/lec10.pdf (presentation)
 	* https://www2.cs.duke.edu/courses/spring07/cps296.2/papers/clarkson-shor.pdf
 * Verification: https://www.spoj.com/problems/CH3D/
 	* https://code.google.com/codejam/contest/6314486/dashboard#s=p3
 */

// using T = ll;
bool above(const P3&a,const P3&b,const P3&c,const P3&p) { 
	return DC(a,b,c,p) > 0; } // is p strictly above plane
void prep(vP3& p) { // rearrange points such that
	shuffle(all(p),rng); // first four are not coplanar
	int dim = 1; 
	FOR(i,1,sz(p)) 
		if (dim == 1) {
			if (p[0] != p[i]) swap(p[1],p[i]), ++dim;
		} else if (dim == 2) {
			if (!collinear(p[0],p[1],p[i])) 
				swap(p[2],p[i]), ++dim;
		} else if (dim == 3) {
			if (!coplanar(p[0],p[1],p[2],p[i]))
				swap(p[3],p[i]), ++dim;
		}
	assert(dim == 4);
}

using F = AR<int,3>; // face
V<F> hull3d(vP3& p) {
	// s.t. first four points form tetra
	prep(p); int N = sz(p); V<F> hull; // triangle for each face
	auto ad = [&](int a, int b, int c) { hull.pb({a,b,c}); }; 
	// +new face to hull
	ad(0,1,2), ad(0,2,1); // initialize hull as first 3 points
	V<vb> in(N,vb(N)); // is zero before each iteration
	FOR(i,3,N) { // incremental construction
		V<F> def, HULL; swap(hull,HULL); 
		// HULL now contains old hull
		auto ins = [&](int a, int b, int c) {
			if (in[b][a]) in[b][a] = 0; // kill reverse face
			else in[a][b] = 1, ad(a,b,c);
		};
		each(f,HULL) {
			if (above(p[f[0]],p[f[1]],p[f[2]],p[i])) 
				F0R(j,3) ins(f[j],f[(j+1)%3],i); 
				// recalc all faces s.t. point is above face
			else def.pb(f); 
		}
		each(t,hull) if (in[t[0]][t[1]]) // edge exposed, 
			in[t[0]][t[1]] = 0, def.pb(t); // add a new face
		swap(hull,def);
	}
	return hull;
}
V<F> hull3dFast(vP3& p) {
	prep(p); int N = sz(p); V<F> hull; 
	vb active; // whether face is active
	V<vi> rvis; // points visible from each face
	V<AR<pi,3>> other; // other face adjacent to each edge of face
	V<vi> vis(N); // faces visible from each point
	auto ad = [&](int a, int b, int c) { 
		hull.pb({a,b,c}); active.pb(1); rvis.eb(); other.eb(); };
	auto ae = [&](int a, int b) { vis[b].pb(a), rvis[a].pb(b); };
	auto abv = [&](int a, int b) {
		F f=hull[a]; return above(p[f[0]],p[f[1]],p[f[2]],p[b]);};
	auto edge = [&](pi e) -> pi { 
		return {hull[e.f][e.s],hull[e.f][(e.s+1)%3]}; };
	auto glue = [&](pi a, pi b) { // link two faces by an edge
		pi x = edge(a); assert(edge(b) == mp(x.s,x.f));
		other[a.f][a.s] = b, other[b.f][b.s] = a;
	}; // ensure face 0 is removed when i=3
	ad(0,1,2), ad(0,2,1); if (abv(1,3)) swap(p[1],p[2]); 
	F0R(i,3) glue({0,i},{1,2-i});
	FOR(i,3,N) ae(abv(1,i),i); // coplanar points go in rvis[0]
	vi label(N,-1);
	FOR(i,3,N) { // incremental construction
		vi rem; each(t,vis[i]) if (active[t]) active[t]=0, rem.pb(t);
		if (!sz(rem)) continue; // hull unchanged
		int st = -1; 
		each(r,rem) F0R(j,3) {
			int o = other[r][j].f;
			if (active[o]) { // create new face!
				int a,b; tie(a,b) = edge({r,j}); ad(a,b,i); st = a;
				int cur = sz(rvis)-1; label[a] = cur; 
				vi tmp; set_union(all(rvis[r]),all(rvis[o]),
									back_inserter(tmp)); 
				// merge sorted vectors ignoring duplicates
				each(x,tmp) if (abv(cur,x)) ae(cur,x);
				/// if no rounding errors then guaranteed that only x>i matters
				glue({cur,0},other[r][j]); // glue old w/ new face
			}
		}
		for (int x = st, y; ; x = y) { // glue new faces together
			int X = label[x]; glue({X,1},{label[y=hull[X][1]],2});
			if (y == st) break;
		}
	}
	V<F> ans; F0R(i,sz(hull)) if (active[i]) ans.pb(hull[i]);
	return ans;
}

/**
 * Description: surface area and volume of polyhedron,
 	* normals to faces must point outwards
 * Source: KACTL
 * Verification: See Hull3D
 */

// #include "Hull3D.h"

pair<T,T> SaVol(vP3 p, V<F> faces) {
	T s = 0, v = 0; 
	each(i,faces) {
		P3 a = p[i[0]], b = p[i[1]], c = p[i[2]];
		s += abs(cross(a,b,c)); v += dot(cross(a,b),c);
	}
	return {s/2,v/6};
}

int main() {
	// read read read
	setIO();
	vP3 v;
	v.pb({-1,1,1});
	v.pb({-1,-1,1});
	v.pb({1,0,-1});
	v.pb({1,1,1});
	v.pb({0,1,-1});
	auto faces = hull3d(v);
	ps(SaVol(v, faces));
	// you should actually read the stuff at the bottom
}

/* stuff you should look for
	* int overflow, array bounds
	* special cases (n=1?)
	* do smth instead of nothing and stay organized
	* WRITE STUFF DOWN
	* DON'T GET STUCK ON ONE APPROACH
*/
