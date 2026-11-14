// Build: g++ -std=c++23 -O3 -march=native -fopenmp -I<damm>/inc bench_ez_update.cc -o bench_ez_update
// Run:   OMP_NUM_THREADS=1 ./bench_ez_update 96 30
// Ez interior update: branchy baseline vs tight row loop vs DAMM composition.
#include <damm.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <cmath>
using T = double;
static T*** alloc3(int a,int b,int c){ // contiguous, rows padded to 8 doubles, 64B aligned
  int cp=c; T* d=(T*)aligned_alloc(64,sizeof(T)*a*b*cp); for(size_t q=0;q<(size_t)a*b*cp;++q) d[q]=std::sin(0.001*q);
  T*** p=new T**[a]; for(int i=0;i<a;++i){p[i]=new T*[b]; for(int j=0;j<b;++j) p[i][j]=d+((size_t)i*b+j)*cp;} return p;}
int main(int argc,char**argv){
  int n=argc>1?atoi(argv[1]):96, steps=argc>2?atoi(argv[2]):50;
  T ***Ez=alloc3(n+1,n+1,n+1),***Hx=alloc3(n+1,n+1,n+1),***Hy=alloc3(n+1,n+1,n+1);
  T ***E2=alloc3(n+1,n+1,n+1),***E3=alloc3(n+1,n+1,n+1);
  std::vector<unsigned char> mat((size_t)(n+1)*(n+1)*(n+1),0);
  T ce1[2]={0.999,0.98}, ce2[2]={0.5,0.25}; T idx=1/1e-9, idy=1/1e-9; idx*=1e-9; idy*=1e-9;
  auto now=[]{return std::chrono::steady_clock::now();};
  auto ms=[](auto a,auto b){return std::chrono::duration<double,std::milli>(b-a).count();};
  // baseline: per-cell region test + material lookup
  int kmin=0,kmax=n+1;
  auto t0=now();
  for(int s=0;s<steps;++s) for(int i=1;i<=n;++i) for(int j=1;j<=n;++j) for(int k=0;k<=n;++k){
    int W = (k>=kmin && k<kmax)?1:(k<kmin?0:2);
    T dy=(Hy[i][j][k]-Hy[i-1][j][k])*idx, dx=(Hx[i][j][k]-Hx[i][j-1][k])*idy, r;
    int m; if(W==1){ r=dy-dx; m=mat[((size_t)i*(n+1)+j)*(n+1)+k]; } else { r=dy-dx; m=0; }
    Ez[i][j][k]*=ce1[m]; Ez[i][j][k]+=ce2[m]*r; }
  auto t1=now();
  // tight: one fused pass per row, scalar coefficients for a uniform run
  for(int s=0;s<steps;++s) for(int i=1;i<=n;++i) for(int j=1;j<=n;++j){
    T* __restrict e=E2[i][j]; const T* a=Hy[i][j]; const T* b=Hy[i-1][j]; const T* c=Hx[i][j]; const T* d=Hx[i][j-1];
    const T c1=ce1[0], c2=ce2[0];
    #pragma omp simd
    for(int k=0;k<=n;++k){ T r=(a[k]-b[k])*idx-(c[k]-d[k])*idy; T v=e[k]*c1; e[k]=v+c2*r; } }
  auto t2=now();
  // DAMM: five whole-plane passes with two scratch planes
  T **t1p=new T*[n], **t2p=new T*[n]; int cp=n+1;
  T* sa=(T*)aligned_alloc(64,sizeof(T)*n*cp); T* sb=(T*)aligned_alloc(64,sizeof(T)*n*cp);
  for(int j=0;j<n;++j){t1p[j]=sa+(size_t)j*cp; t2p[j]=sb+(size_t)j*cp;}
  using namespace damm;
  for(int s=0;s<steps;++s) for(int i=1;i<=n;++i){
    size_t M=n, N=n+1;
    scalar::fused_union<FusionPolicy::UNION_FIRST,T,std::minus<>,std::multiplies<>>(Hy[i]+1,Hy[i-1]+1,idx,t1p,M,N);
    scalar::fused_union<FusionPolicy::UNION_FIRST,T,std::minus<>,std::multiplies<>>(Hx[i]+1,Hx[i],idy,t2p,M,N);
    matrix::unite<T,std::minus<>>(t1p,t2p,t1p,M,N);
    scalar::unite<T,std::multiplies<>>(E3[i]+1,ce1[0],E3[i]+1,M,N);
    scalar::fused_union<FusionPolicy::FUSION_FIRST,T,std::plus<>,std::multiplies<>>(E3[i]+1,ce2[0],t1p,E3[i]+1,M,N);
  }
  auto t3=now();
  double e12=0,e13=0; for(int i=1;i<=n;++i)for(int j=1;j<=n;++j)for(int k=0;k<=n;++k){e12=std::fmax(e12,std::fabs(Ez[i][j][k]-E2[i][j][k])); e13=std::fmax(e13,std::fabs(Ez[i][j][k]-E3[i][j][k]));}
  double cells=(double)n*n*(n+1)*steps;
  printf("n=%d steps=%d\n baseline %.1f ms (%.2f ns/cell)\n tight    %.1f ms (%.2f ns/cell)  maxdiff %.3g\n damm     %.1f ms (%.2f ns/cell)  maxdiff %.3g\n",
    n,steps,ms(t0,t1),ms(t0,t1)*1e6/cells,ms(t1,t2),ms(t1,t2)*1e6/cells,e12,ms(t2,t3),ms(t2,t3)*1e6/cells,e13);
}
