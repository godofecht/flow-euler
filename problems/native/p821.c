#include <stdint.h>
#include <stdlib.h>
typedef long long i64;
static i64 coprime6(i64 m){ if(m<=0) return 0; return m - m/2 - m/3 + m/6; }
static int cmp_i64(const void*a,const void*b){ i64 x=*(const i64*)a,y=*(const i64*)b; return (x>y)-(x<y); }
static int best_covered_small(i64 limit, i64 *vals, int m){
    int best=0;
    for(int mask=0; mask<(1<<m); mask++){
        char usedS[64]={0}, usedD[64]={0}, usedT[64]={0};
        /* map values to indices 0..m-1 for membership among smooth; also track 2v,3v if in set */
        int ok=1; int covered=0;
        i64 S[32],D[32],T[32]; int ns=0,nd=0,nt=0;
        for(int i=0;i<m;i++) if(mask>>i & 1){
            i64 v=vals[i], dv=2*v, tv=3*v;
            for(int j=0;j<nd;j++) if(S[j]==v||D[j]==v||T[j]==v) {ok=0;break;}
            /* check disjoint */
            for(int j=0;j<ns;j++) if(S[j]==dv||S[j]==tv) ok=0;
            for(int j=0;j<nd;j++) if(D[j]==v||D[j]==tv) ok=0;
            for(int j=0;j<nt;j++) if(T[j]==v||T[j]==dv) ok=0;
            if(!ok) break;
            S[ns++]=v; D[nd++]=dv; T[nt++]=tv;
        }
        if(!ok) continue;
        /* count unique <= limit */
        i64 u[96]; int nu=0;
        for(int i=0;i<ns;i++) u[nu++]=S[i];
        for(int i=0;i<nd;i++) u[nu++]=D[i];
        for(int i=0;i<nt;i++) u[nu++]=T[i];
        qsort(u,(size_t)nu,sizeof(i64),cmp_i64);
        for(int i=0;i<nu;i++) if(u[i]<=limit && (i==0||u[i]!=u[i-1])) covered++;
        if(covered>best) best=covered;
    }
    return best;
}
long long pe_solve(void){
    i64 n=10000000000000000LL;
    i64 smooth[4096]; int ns=0;
    for(i64 p2=1; p2<=n; p2*=2){
        for(i64 p3=1; p2*p3<=n; p3*=3) smooth[ns++]=p2*p3;
        if(p2 > n/2) break;
    }
    qsort(smooth,(size_t)ns,sizeof(i64),cmp_i64);
    /* unique */
    int w=0; for(int i=0;i<ns;i++) if(i==0||smooth[i]!=smooth[i-1]) smooth[w++]=smooth[i]; ns=w;
    int small_H[64]={0};
    for(int i=0;i<ns && smooth[i]<=48;i++){
        i64 vals[32]; int m=0;
        for(int j=0;j<=i;j++) vals[m++]=smooth[j];
        small_H[i]=best_covered_small(smooth[i], vals, m);
    }
    i64 holes[256]; int nh=0;
    for(i64 x=6; ; ){ /* fixed list */
        i64 init[]={6,24,54}; int done=1;
        for(int i=0;i<3;i++) if(init[i]<=n){ holes[nh++]=init[i]; }
        break;
    }
    for(i64 x=384; x<=n; x*=8) holes[nh++]=x;
    for(i64 x=243; x<=n; x*=27) holes[nh++]=x;
    qsort(holes,(size_t)nh,sizeof(i64),cmp_i64);
    int hole_prefix[4096]; int j=0,cnt=0;
    for(int i=0;i<ns;i++){
        while(j<nh && holes[j]<=smooth[i]){ cnt++; j++; }
        hole_prefix[i]=cnt;
    }
    i64 ans=0;
    for(int i=0;i<ns;i++){
        i64 L=smooth[i];
        i64 R = (i+1==ns)? n : (smooth[i+1]-1 < n ? smooth[i+1]-1 : n);
        i64 k_low = n/(R+1)+1;
        i64 k_high = n/L;
        if(k_low>k_high) continue;
        i64 count_k = coprime6(k_high)-coprime6(k_low-1);
        i64 H;
        if(L<=48) H=small_H[i];
        else H=(i+1)-hole_prefix[i];
        ans += H * count_k;
    }
    return ans;
}
