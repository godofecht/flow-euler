#include <stdint.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
enum { MOD = 1000000007LL };
static i64 modpow(i64 b,i64 e){i64 r=1;b%=MOD;while(e){if(e&1)r=r*b%MOD;b=b*b%MOD;e>>=1;}return r;}
long long pe_solve(void){
    int n=1000000, k=999983;
    int B=0; for(int t=n;t;t>>=1) B++;
    int size=1<<B;
    unsigned char *is_prime=malloc((size_t)n+1); memset(is_prime,1,(size_t)n+1);
    is_prime[0]=is_prime[1]=0;
    for(int p=2;(i64)p*p<=n;p++) if(is_prime[p]) for(int j=p*p;j<=n;j+=p) is_prime[j]=0;
    i64 *a=calloc((size_t)size,sizeof(i64));
    i64 *s=calloc((size_t)size,sizeof(i64));
    for(int p=2;p<=n;p++) if(is_prime[p]){
        a[p]=1; s[p]=(__builtin_popcount((unsigned)p)&1)?-1:1;
    }
    for(int i=0;i<B;i++){
        int step=1<<i, jump=step<<1;
        for(int base=0;base<size;base+=jump)
            for(int m=base+step;m<base+jump;m++) a[m]+=a[m-step];
    }
    for(int i=0;i<B;i++){
        int step=1<<i, jump=step<<1;
        for(int base=0;base<size;base+=jump)
            for(int m=base;m<base+step;m++) s[m]+=s[m+step];
    }
    i64 res=0;
    for(int m=0;m<=n;m++){
        if(!a[m]) continue;
        i64 val=modpow(a[m],k);
        i64 coeff=s[m];
        if(__builtin_popcount((unsigned)m)&1) coeff=-coeff;
        res=(res+val*coeff)%MOD;
    }
    if(res<0) res+=MOD;
    free(is_prime); free(a); free(s);
    return res;
}
