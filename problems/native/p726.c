#include <stdint.h>
#include <stdlib.h>
typedef long long i64;
enum { MOD = 1000000033 };
static i64 modpow(i64 b,i64 e){i64 r=1;b%=MOD;while(e){if(e&1)r=r*b%MOD;b=b*b%MOD;e>>=1;}return r;}
long long pe_solve(void){
    int n=10000;
    i64 *inv_odd=malloc((size_t)(n+1)*sizeof(i64));
    for(int k=1;k<=n;k++) inv_odd[k]=modpow(2*k-1, MOD-2);
    i64 cur_f=1, total=1, curN=1, pow2=2, mers_prefix=1, odd_inv_prefix=1;
    for(int layer=2; layer<=n; layer++){
        i64 start=curN+1, end=curN+layer;
        for(i64 x=start;x<=end;x++) cur_f=cur_f*x%MOD;
        curN=end;
        pow2=pow2*2%MOD;
        mers_prefix=mers_prefix*((pow2-1+MOD)%MOD)%MOD;
        odd_inv_prefix=odd_inv_prefix*inv_odd[layer]%MOD;
        cur_f=cur_f*mers_prefix%MOD*odd_inv_prefix%MOD;
        total=(total+cur_f)%MOD;
    }
    free(inv_odd);
    return total;
}
