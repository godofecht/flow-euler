#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
enum { MOD = 1000000007LL };
long long pe_solve(void){
    i64 k=100000000LL;
    i64 mod=MOD;
    /* lcm(1..k) mod = product p^{floor(log_p k)} */
    /* = 2 * product odd primes p * p^{e-1} for e=floor(log_p k)>=1 ... */
    /* Port of lcm_1_to_k_mod from python */
    i64 res=1;
    /* handle 2 */
    i64 p2=2;
    while(p2*2<=k) p2*=2;
    res=p2%mod;
    i64 sqrt_k=(i64)sqrt((double)k);
    while((sqrt_k+1)*(sqrt_k+1)<=k)sqrt_k++;
    /* sieve small primes */
    char *comp=calloc((size_t)sqrt_k+1,1);
    for(i64 i=2;i*i<=sqrt_k;i++) if(!comp[i]) for(i64 j=i*i;j<=sqrt_k;j+=i) comp[j]=1;
    /* for each small prime, mark in segments and multiply */
    /* Simpler approach: sieve all primes up to k with segmented sieve */
    i64 low=3;
    i64 SEG=1000000;
    char *seg=malloc((size_t)SEG);
    while(low<=k){
        i64 high=low+SEG; if(high>k+1) high=k+1;
        memset(seg,0,(size_t)(high-low));
        for(i64 p=3;p<=sqrt_k;p++){
            if(comp[p]) continue;
            i64 start=((low+p-1)/p)*p;
            if(start<p*p) start=p*p;
            if((start&1)==0) start+=p; /* keep odd */
            for(i64 m=start;m<high;m+=2*p) seg[m-low]=1;
        }
        for(i64 n=low;n<high;n+=2){
            if(!seg[n-low]){
                res=res*n%mod;
                if(n<=sqrt_k){
                    i64 power=n*n;
                    while(power<=k){ res=res*n%mod; power*=n; }
                }
            }
        }
        low=high;
        if((low&1)==0) low++;
    }
    free(comp); free(seg);
    return (2*res)%mod; /* f(k)=2*lcm */
}
