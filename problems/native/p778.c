#include <stdint.h>
#include <string.h>
typedef long long i64;
enum { MOD = 1000000009LL };
static void mat_mul(i64 C[10][10], i64 A[10][10], i64 B[10][10]){
    i64 T[10][10]; memset(T,0,sizeof T);
    for(int i=0;i<10;i++) for(int k=0;k<10;k++) if(A[i][k])
        for(int j=0;j<10;j++) T[i][j]=(T[i][j]+A[i][k]*B[k][j])%MOD;
    memcpy(C,T,sizeof T);
}
static void mat_pow(i64 R[10][10], i64 A0[10][10], i64 exp){
    i64 A[10][10]; memcpy(A,A0,sizeof A);
    memset(R,0,sizeof(i64)*100); for(int i=0;i<10;i++) R[i][i]=1;
    while(exp){ if(exp&1) mat_mul(R,R,A); mat_mul(A,A,A); exp>>=1; }
}
static void digit_counts(i64 n, int pos, i64 counts[10]){
    i64 base=1; for(int i=0;i<pos;i++) base*=10;
    i64 higher=n/(base*10), cur=(n/base)%10, lower=n%base;
    for(int d=0;d<10;d++) counts[d]=higher*base;
    for(int d=0;d<cur;d++) counts[d]+=base;
    counts[cur]+=lower+1;
}
long long pe_solve(void){
    i64 R=234567, M=765432;
    int max_digits=0; for(i64 t=M;t;t/=10) max_digits++;
    i64 ans=0, pow10=1;
    for(int pos=0;pos<max_digits;pos++){
        i64 counts[10]; digit_counts(M,pos,counts);
        i64 A[10][10]; memset(A,0,sizeof A);
        for(int s=0;s<10;s++) for(int d=0;d<10;d++){
            int t=(s*d)%10; A[s][t]=(A[s][t]+counts[d])%MOD;
        }
        i64 P[10][10]; mat_pow(P,A,R);
        i64 digit_sum=0;
        for(int digit=0;digit<10;digit++) digit_sum=(digit_sum+digit*P[1][digit])%MOD;
        ans=(ans+digit_sum*pow10)%MOD;
        pow10=pow10*10%MOD;
    }
    return ans;
}
