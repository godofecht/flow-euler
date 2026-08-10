#include <stdint.h>
#include <stdlib.h>
typedef long long i64;
typedef __int128 i128;
static i64 pow_i(i64 p,int e){ i64 r=1; for(int i=0;i<e;i++) r*=p; return r; }
static i64 a_prime_power(i64 p,int e){
    if(e==0) return 1;
    if(e&1){ int k=e>>1; return pow_i(p,e-1)*(2*pow_i(p,k+1)-1); }
    int k=e>>1; return pow_i(p,e-1)*((p+1)*pow_i(p,k)-1);
}
long long pe_solve(void){
    i64 N=12345678;
    i64 limit=N/2;
    i64 *spf=calloc((size_t)(limit+1),sizeof(i64));
    i64 *primes=malloc((size_t)(limit+1)*sizeof(i64)); i64 pc=0;
    for(i64 i=2;i<=limit;i++){
        if(!spf[i]){ spf[i]=i; primes[pc++]=i; }
        for(i64 j=0;j<pc;j++){
            i64 p=primes[j]; i64 x=p*i; if(x>limit||p>spf[i]) break; spf[x]=p; if(p==spf[i]) break;
        }
    }
    int max_a=64; i64 c2[65]={0};
    for(int a=1;a<=64;a++){
        if((1LL<<a)>N) break;
        c2[a]=a_prime_power(2,a)-(1LL<<a);
    }
    i64 odd_cnt=(N+1)/2;
    i64 total=-(odd_cnt*odd_cnt);
    for(i64 m=1;m<=limit;m+=2){
        i64 a_m=1;
        if(m!=1){
            i64 n=m; a_m=1;
            while(n>1){
                i64 p=spf[n]; int e=0; while(n%p==0){ n/=p; e++; }
                a_m*=a_prime_power(p,e);
            }
        }
        int a=1; i64 n=m<<1;
        while(n<=N){ total+=a_m*c2[a]; a++; n<<=1; }
    }
    free(spf); free(primes);
    return total;
}
