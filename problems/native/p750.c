#include <stdint.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
enum { INF = (i64)4e18 };
long long pe_solve(void){
    int N=976;
    int mod=N+1;
    int *pos=calloc((size_t)(N+1),sizeof(int));
    int x=1;
    for(int i=1;i<=N;i++){
        x=(x*3)%mod;
        if(x==0||x>N||pos[x]) { free(pos); return -1; }
        pos[x]=i;
    }
    i64 *dp=malloc((size_t)(N+1)*(N+1)*sizeof(i64));
    #define DP(l,r) dp[(l)*(N+1)+(r)]
    for(int i=0;i<=N;i++) for(int j=0;j<=N;j++) DP(i,j)=0;
    for(int r=2;r<=N;r++){
        for(int l=r-1;l>=1;l--){
            i64 best=INF;
            for(int k=l;k<r;k++){
                i64 cost=DP(l,k)+DP(k+1,r)+llabs((i64)pos[k]-pos[r]);
                if(cost<best) best=cost;
            }
            DP(l,r)=best;
        }
    }
    i64 ans=DP(1,N);
    free(dp); free(pos);
    return ans;
}
