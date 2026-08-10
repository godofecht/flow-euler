#include <math.h>
#include <stdio.h>
#include <stdint.h>
typedef long long i64;
typedef __int128 i128;
static i64 modpow(i64 b, i64 e, i64 m){
    i64 r=1%m; b%=m;
    while(e>0){ if(e&1) r=(i64)((i128)r*b%m); b=(i64)((i128)b*b%m); e>>=1; }
    return r;
}
void pe_solve_print(void){
    i64 n=10000000000000000LL;
    int limit=(int)(log((double)n)/log(3.0));
    long double total=0;
    for(int i=1;i<=limit;i++){
        i64 div=1; for(int j=0;j<i;j++) div*=3;
        i64 exp=n-div-1;
        i64 temp=modpow(10,exp,div);
        total += (long double)temp / (long double)div;
    }
    long double frac = total - floorl(total);
    /* Match Python: str(frac)[2:12] — truncate, don't round. */
    unsigned long long digits = (unsigned long long)(frac * 10000000000.0L);
    printf("%010llu\n", digits);
}
