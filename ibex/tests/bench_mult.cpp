// Tiempo por producto de AffineMain<AF_fAF2> con af2, raa (k=0), sh y exact, segun n.
// Factores con error (lin + 0.3*lin*lin) o sin error (lin). [parche 0016]
#include <bits/stdc++.h>
#include "ibex.h"
#include "ibex_Affine.h"
#include "ibex_KNomial.h"
#include "ibex_ExactMult.h"
using namespace std; using namespace ibex;
typedef AffineMain<AF_fAF2> AF;
static mt19937_64 rng(3);
static double U(double a,double b){return uniform_real_distribution<double>(a,b)(rng);}
static AF lin(const vector<AF>& X){ AF r(Interval(U(-2,2))); for(size_t i=0;i<X.size();i++){ AF t=X[i]; t*=U(-1,1); r+=t; } return r; }
int main(){
  const AffineMode modes[4]={AFFINE_AF2,AFFINE_RAA,AFFINE_SH,AFFINE_EXACT}; const char* names[4]={"af2","raa","sh","exact"};
  printf("n,errores,af2_us,raa_us,sh_us,exact_us,exact/sh,exact/af2\n");
  for(int n: {2,5,10,20,50,100}) for(int witherr=0; witherr<2; witherr++){
    IntervalVector Bx(n, Interval(-1,1)); AffineVarMainVector<AF_fAF2> V(Bx); vector<AF> X; for(int i=0;i<n;i++) X.push_back(AF(V[i]));
    const int P=64; vector<AF> A,B;
    affine_mode()=AFFINE_AF2;
    for(int p=0;p<P;p++){ AF a=lin(X), b=lin(X); if(witherr){ AF t=lin(X)*lin(X); t*=0.3; a+=t; AF s=lin(X)*lin(X); s*=0.3; b+=s; } A.push_back(a); B.push_back(b); }
    double t[4];
    for(int m=0;m<4;m++){
      affine_mode()=modes[m];
      const int R = max(20, 40000/n);
      double sink=0; auto t0=chrono::steady_clock::now();
      for(int r=0;r<R;r++) for(int p=0;p<P;p++){ AF c=A[p]*B[p]; sink+=c.err(); }
      auto t1=chrono::steady_clock::now();
      t[m]=chrono::duration<double,micro>(t1-t0).count()/(double(R)*P);
      if(sink<0) printf("!");
    }
    affine_mode()=AFFINE_AF2;
    printf("%d,%s,%.3f,%.3f,%.3f,%.3f,%.2f,%.2f\n",n,witherr?"si":"no",t[0],t[1],t[2],t[3],t[3]/t[2],t[3]/t[0]);
  }
}
