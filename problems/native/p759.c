#include <stdint.h>
#include <string.h>
typedef long long i64;
enum { MOD = 1000000007LL, MAXB = 60 };
typedef i64 Mat[3][3];
static void zero(Mat a){ memset(a,0,sizeof(Mat)); }
static void addm(Mat o, Mat a, Mat b){
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) o[i][j]=(a[i][j]+b[i][j])%MOD;
}
static void shift_range(Mat out, Mat mat, i64 p){
    p%=MOD; i64 p2=p*p%MOD, two_p=2*p%MOD;
    i64 mats[3][3];
    for(int t=0;t<3;t++){
        i64 s0=mat[t][0]%MOD,s1=mat[t][1]%MOD,s2=mat[t][2]%MOD;
        mats[t][0]=s0; mats[t][1]=(p*s0+s1)%MOD; mats[t][2]=(p2*s0+two_p*s1+s2)%MOD;
    }
    int coeffs[3][3]={{1,0,0},{1,1,0},{1,2,1}};
    zero(out);
    for(int j=0;j<3;j++) for(int d=0;d<3;d++)
        out[j][d]=(coeffs[j][0]*mats[0][d]+coeffs[j][1]*mats[1][d]+coeffs[j][2]*mats[2][d])%MOD;
}
static Mat full[MAXB+1];
static void precompute(int max_bits){
    zero(full[0]); full[0][0][0]=1;
    for(int m=1;m<=max_bits;m++){
        Mat sh; shift_range(sh, full[m-1], 1LL<<(m-1));
        addm(full[m], full[m-1], sh);
    }
}
static void calc_upto(Mat out, i64 n){
    if(n<0){ zero(out); return; }
    if(n==0){ memcpy(out,full[0],sizeof(Mat)); return; }
    int k=63-__builtin_clzll((unsigned long long)n);
    i64 p=1LL<<k;
    if(n==p-1){ memcpy(out,full[k],sizeof(Mat)); return; }
    Mat left, right_in, sh;
    memcpy(left, full[k], sizeof(Mat));
    calc_upto(right_in, n-p);
    shift_range(sh, right_in, p);
    addm(out, left, sh);
}
long long pe_solve(void){
    i64 n=10000000000000000LL;
    precompute(60);
    Mat m; calc_upto(m,n);
    return m[2][2]%MOD;
}
