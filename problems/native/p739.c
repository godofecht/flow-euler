#include <stdint.h>
#include <stdlib.h>
typedef long long i64;
enum { MOD = 1000000007LL, INV2 = (MOD + 1) / 2 };
static i64 modpow(i64 b, i64 e){i64 r=1;b%=MOD;while(e){if(e&1)r=r*b%MOD;b=b*b%MOD;e>>=1;}return r;}
static void inverses_consecutive(i64 start, int length, i64 *invs){
    i64 *pref = malloc((size_t)(length+1)*sizeof(i64));
    pref[0]=1; i64 x=start;
    for(int i=0;i<length;i++){ pref[i+1]=pref[i]*x%MOD; x++; }
    i64 inv_total = modpow(pref[length], MOD-2);
    x = start+length-1;
    for(int i=length-1;i>=0;i--){
        invs[i] = inv_total * pref[i] % MOD * INV2 % MOD;
        inv_total = inv_total * x % MOD;
        x--;
    }
    free(pref);
}
long long pe_solve(void){
    i64 n=100000000LL;
    i64 m=n-1;
    i64 a0=1,a1=3,a2=7,a3=21;
    i64 steps=m-3;
    i64 c0=2,c1=26,c2=62,c3=50;
    i64 denom=4;
    const int BLOCK=200000;
    i64 *invs=malloc((size_t)BLOCK*sizeof(i64));
    while(steps>0){
        int L = steps>BLOCK?BLOCK:(int)steps;
        inverses_consecutive(denom, L, invs);
        for(int j=0;j<L;j++){
            i64 t=(c3*a3%MOD - c0*a0%MOD - c1*a1%MOD - c2*a2%MOD)%MOD;
            if(t<0)t+=MOD;
            i64 a4=t*invs[j]%MOD;
            a0=a1;a1=a2;a2=a3;a3=a4;
            c0+=4;c1+=23;c2+=22;c3+=15; denom++;
        }
        steps-=L;
    }
    free(invs);
    return a3;
}
