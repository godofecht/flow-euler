#include <stdint.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
enum { MOD = 1001001011 };
long long pe_solve(void){
    int n=20; int N=1<<n;
    int *succ=malloc((size_t)N*sizeof(int));
    int *indeg=calloc((size_t)N,sizeof(int));
    int mask=(1<<(n-1))-1, shift=n-3;
    for(int s=0;s<N;s++){
        int t=s>>shift;
        int newbit=((t>>2)&1)&(((t>>1)&1)^(t&1));
        int ns=((s&mask)<<1)|newbit;
        succ[s]=ns; indeg[ns]++;
    }
    i64 *acc0=malloc((size_t)N*sizeof(i64));
    i64 *acc1=malloc((size_t)N*sizeof(i64));
    for(int i=0;i<N;i++){ acc0[i]=1; acc1[i]=1; }
    unsigned char *in_cycle=malloc((size_t)N); memset(in_cycle,1,(size_t)N);
    int *q=malloc((size_t)N*sizeof(int)); int qh=0,qt=0;
    for(int i=0;i<N;i++) if(indeg[i]==0) q[qt++]=i;
    while(qh<qt){
        int u=q[qh++];
        if(!in_cycle[u]) continue;
        in_cycle[u]=0;
        int p=succ[u];
        i64 dp0=acc0[u], dp1=acc1[u];
        acc0[p]=acc0[p]*((dp0+dp1)%MOD)%MOD;
        acc1[p]=acc1[p]*dp0%MOD;
        if(--indeg[p]==0) q[qt++]=p;
    }
    unsigned char *visited=calloc((size_t)N,1);
    i64 ans=1;
    int *cycle=malloc((size_t)N*sizeof(int));
    for(int v=0;v<N;v++) if(in_cycle[v] && !visited[v]){
        int k=0,u=v;
        while(!visited[u]){ visited[u]=1; cycle[k++]=u; u=succ[u]; }
        i64 case1, case2;
        {
            i64 prev0=acc0[cycle[0]]%MOD, prev1=0;
            for(int i=1;i<k;i++){
                i64 cur0=((prev0+prev1)%MOD)*(acc0[cycle[i]]%MOD)%MOD;
                i64 cur1=prev0*(acc1[cycle[i]]%MOD)%MOD;
                prev0=cur0; prev1=cur1;
            }
            case1=(prev0+prev1)%MOD;
        }
        {
            i64 prev0=0, prev1=acc1[cycle[0]]%MOD;
            for(int i=1;i<k;i++){
                i64 cur0=((prev0+prev1)%MOD)*(acc0[cycle[i]]%MOD)%MOD;
                i64 cur1=prev0*(acc1[cycle[i]]%MOD)%MOD;
                prev0=cur0; prev1=cur1;
            }
            case2=prev0;
        }
        ans=ans*((case1+case2)%MOD)%MOD;
    }
    free(succ);free(indeg);free(acc0);free(acc1);free(in_cycle);free(q);free(visited);free(cycle);
    return ans;
}
