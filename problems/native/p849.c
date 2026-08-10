#include <stdint.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
enum { MOD = 1000000007LL };
long long pe_solve(void){
    int n=100;
    int vmax=2*(n-1), off=vmax, V=2*vmax+1;
    int *maxh=malloc((size_t)(n+1)*sizeof(int));
    for(int i=0;i<=n;i++) maxh[i]=2*i*(n-i);
    /* dp as flat: (mh+1)*V */
    int mh=maxh[1];
    i64 *dp=calloc((size_t)(mh+1)*V, sizeof(i64));
    for(int d=0;d<=vmax;d++) dp[d*V + (d+off)]=1;
    for(int step=1; step<n; step++){
        int mh_cur=maxh[step], mh_next=maxh[step+1];
        i64 *nxt=calloc((size_t)(mh_next+1)*V, sizeof(i64));
        for(int h=0;h<=mh_cur;h++){
            i64 *row=dp+h*V;
            i64 cum=0;
            int base=h-off-4;
            int vi_start=-base; if(vi_start<4) vi_start=4; if(vi_start>V) vi_start=V;
            int vi_end=mh_next-base; if(vi_end>V-1) vi_end=V-1;
            for(int vi=0; vi<vi_start; vi++){ cum+=row[vi]; if(cum>=MOD) cum-=MOD; }
            if(vi_start<=vi_end){
                for(int vi=vi_start; vi<=vi_end; vi++){
                    cum+=row[vi]; if(cum>=MOD) cum-=MOD;
                    int h2=base+vi; int idx=vi-4;
                    i64 val=nxt[h2*V+idx]+cum; if(val>=MOD) val-=MOD; nxt[h2*V+idx]=val;
                }
                for(int vi=vi_end+1; vi<V; vi++){ cum+=row[vi]; if(cum>=MOD) cum-=MOD; }
            } else {
                for(int vi=vi_start; vi<V; vi++){ cum+=row[vi]; if(cum>=MOD) cum-=MOD; }
            }
            i64 total=cum;
            if(total){
                int base_h=h-off;
                for(int vi=V-4; vi<V; vi++){
                    int h2=base_h+vi;
                    if(0<=h2 && h2<=mh_next){
                        i64 val=nxt[h2*V+vi]+total; if(val>=MOD) val-=MOD; nxt[h2*V+vi]=val;
                    }
                }
            }
        }
        free(dp); dp=nxt;
    }
    i64 ans=0; for(int vi=0; vi<V; vi++){ ans+=dp[0*V+vi]; if(ans>=MOD) ans-=MOD; }
    free(dp); free(maxh);
    return ans;
}
