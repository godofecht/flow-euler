#include <stdint.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
typedef __int128 i128;
typedef struct { i64 k,cap; i128 val; int used; } Ent;
enum { HSIZE = 1<<22 };
static Ent *H;
static i64 *F;
static i64 sum_all(i64 k){ if(k<=0) return 0; return F[k+2]-2; }
static i128 count(i64 k, i64 cap){
    if(cap<0) return 0;
    if(k==0) return 1;
    i64 key = (k<<48) ^ (cap & ((1LL<<48)-1)); /* rough */
    uint64_t h = ((uint64_t)k*11400714819323198485ull) ^ (uint64_t)cap;
    for(int probe=0; probe<64; probe++){
        Ent *e=&H[(h+probe)&(HSIZE-1)];
        if(e->used && e->k==k && e->cap==cap) return e->val;
        if(!e->used){
            i128 res;
            if(cap >= sum_all(k)){
                if(k>=127) res = ((i128)1)<<k; /* won't happen for our n */
                else res = ((i128)1)<<k;
            } else if(cap < F[k]) res = count(k-1,cap);
            else res = count(k-1,cap)+count(k-1,cap-F[k]);
            e->used=1; e->k=k; e->cap=cap; e->val=res;
            return res;
        }
    }
    /* fallback no cache */
    if(cap >= sum_all(k)) return ((i128)1)<<k;
    if(cap < F[k]) return count(k-1,cap);
    return count(k-1,cap)+count(k-1,cap-F[k]);
}
long long pe_solve(void){
    i64 n=10000000000000LL;
    F=calloc(200,sizeof(i64));
    F[1]=1; F[2]=2; int len=2;
    while(F[len]<=n){ len++; F[len]=F[len-1]+F[len-2]; }
    int k=len-1;
    while(len<=k+2){ len++; F[len]=F[len-1]+F[len-2]; }
    H=calloc(HSIZE,sizeof(Ent));
    i128 ans=count(k,n);
    free(H); free(F);
    return (long long)ans;
}
