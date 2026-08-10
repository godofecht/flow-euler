#include <stdint.h>
#include <string.h>
typedef long long i64;
enum { MOD=1000000007LL, N=81 };
static int idx(int a,int b,int c,int cur){ return ((a*3+b)*3+c)*3+cur; }
static void mat_mul(i64 C[N][N], i64 A[N][N], i64 B[N][N]){
    i64 T[N][N]; memset(T,0,sizeof T);
    for(int i=0;i<N;i++) for(int k=0;k<N;k++) if(A[i][k])
        for(int j=0;j<N;j++) if(B[k][j]) T[i][j]=(T[i][j]+A[i][k]*B[k][j])%MOD;
    memcpy(C,T,sizeof T);
}
static void mat_pow(i64 R[N][N], i64 A0[N][N], i64 e){
    i64 A[N][N]; memcpy(A,A0,sizeof A);
    memset(R,0,sizeof(i64)*N*N); for(int i=0;i<N;i++) R[i][i]=1;
    while(e){ if(e&1) mat_mul(R,A,R); e>>=1; if(e) mat_mul(A,A,A); }
}
static void mat_vec(i64 out[N], i64 M[N][N], i64 v[N]){
    for(int i=0;i<N;i++){ i64 s=0; for(int j=0;j<N;j++) if(M[i][j]) s=(s+M[i][j]*v[j])%MOD; out[i]=s; }
}
static void build(i64 M[N][N], i64 w0,i64 w1,i64 w2){
    memset(M,0,sizeof(i64)*N*N);
    i64 w[3]={w0,w1,w2};
    for(int a=0;a<3;a++) for(int b=0;b<3;b++) for(int c=0;c<3;c++) for(int cur=0;cur<3;cur++){
        int frm=idx(a,b,c,cur);
        for(int r=0;r<3;r++){
            int nxt=(cur+r)%3;
            int a2=a,b2=b,c2=c;
            if(nxt==0) a2=(a2+1)%3; else if(nxt==1) b2=(b2+1)%3; else c2=(c2+1)%3;
            int to=idx(a2,b2,c2,nxt);
            M[to][frm]=(M[to][frm]+w[r])%MOD;
        }
    }
}
long long pe_solve(void){
    int d=100000;
    i64 lead[N][N], step[N][N];
    build(lead,3,3,3); build(step,4,3,3);
    int good[N]={0};
    for(int a=0;a<3;a++) for(int b=0;b<3;b++) for(int c=0;c<3;c++){
        int k=(a==2)+(b==2)+(c==2);
        if(k%3==0) for(int cur=0;cur<3;cur++) good[idx(a,b,c,cur)]=1;
    }
    i64 v[N]={0}, v2[N];
    v[idx(1,0,0,0)]=1;
    mat_vec(v2, lead, v); memcpy(v,v2,sizeof v);
    if(d>1){
        i64 P[N][N]; mat_pow(P, step, d-1);
        mat_vec(v2, P, v); memcpy(v,v2,sizeof v);
    }
    i64 ans=0; for(int i=0;i<N;i++) if(good[i]) ans+=v[i];
    return ans%MOD;
}
