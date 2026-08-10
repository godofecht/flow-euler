#include <stdint.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
static int bit_length(i64 x){ int c=0; while(x){x>>=1;c++;} return c; }
static int search_bit_limit(int rank){
    /* heuristic from growth; expand until sieve finds enough */
    for(int bl=20; bl<40; bl++){
        /* rough: number of xor-primes ~ similar to primes; use generous */
        if((1LL<<bl)/bl > (i64)rank*3) return bl;
    }
    return 32;
}
long long pe_solve(void){
    int rank=5000000;
    if(rank==1) return 2;
    int bit_limit;
    /* find sufficient bit_limit by retrying */
    for(bit_limit=24; bit_limit<=34; bit_limit++){
        i64 limit=1LL<<bit_limit;
        unsigned char *mark=calloc((size_t)(limit>>1),1);
        if(!mark) continue;
        mark[0]=1;
        int found=1;
        i64 ans=-1;
        for(i64 base=3; base<limit; base+=2){
            if(mark[base>>1]) continue;
            found++;
            if(found==rank){ ans=base; break; }
            int degree=bit_length(base)-1;
            int max_cofactor_degree=bit_limit-degree-1;
            for(int cofactor_degree=degree; cofactor_degree<=max_cofactor_degree; cofactor_degree++){
                i64 product=(base<<cofactor_degree)^base;
                if(product<limit) mark[product>>1]=1;
                i64 variants=1LL<<(cofactor_degree-1);
                for(i64 n=1;n<variants;n++){
                    int toggled_bit = __builtin_ctzll((unsigned long long)(n & -n)) + 1; /* bit_length of lsb */
                    /* python: (n & -n).bit_length() is position of LSB starting at 1 */
                    product ^= base << toggled_bit;
                    if(product>=0 && product<limit) mark[product>>1]=1;
                }
            }
        }
        free(mark);
        if(ans>=0) return ans;
    }
    return -1;
}
