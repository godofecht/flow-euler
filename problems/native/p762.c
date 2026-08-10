#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
enum { MOD = 1000000000LL };
typedef struct { int a[4]; } State;
static int pop4(int m){ return (m&1)+((m>>1)&1)+((m>>2)&1)+((m>>3)&1); }
void pe_solve_print(void){
    State states[70]; int idx[4][4][4][4]; memset(idx,-1,sizeof idx); int S=0;
    for(int a0=0;a0<4;a0++) for(int a1=0;a1<4-a0;a1++) for(int a2=0;a2<4-a0-a1;a2++) for(int a3=0;a3<4-a0-a1-a2;a3++){
        states[S]=(State){{a0,a1,a2,a3}}; idx[a0][a1][a2][a3]=S++;
    }
    int terminal=idx[0][0][0][0];
    int *to_nt_v[70]; int *to_nt_w[70]; int to_nt_n[70]={0};
    int *to_t_w[70]; int to_t_n[70]={0};
    for(int u=0;u<S;u++){ to_nt_v[u]=malloc(64*sizeof(int)); to_nt_w[u]=malloc(64*sizeof(int)); to_t_w[u]=malloc(16*sizeof(int)); }
    for(int u=0;u<S;u++){
        if(u==terminal) continue;
        int *s=states[u].a;
        int t0=s[0]+s[3], t1=s[1]+s[0], t2=s[2]+s[1], t3=s[3]+s[2];
        for(int mask=0;mask<16;mask++){
            int b0=mask&1,b1=(mask>>1)&1,b2=(mask>>2)&1,b3=(mask>>3)&1;
            int n0=t0-b0,n1=t1-b1,n2=t2-b2,n3=t3-b3;
            if(n0<0||n1<0||n2<0||n3<0||n0>3||n1>3||n2>3||n3>3||n0+n1+n2+n3>3) continue;
            int v=idx[n0][n1][n2][n3]; if(v<0) continue;
            int w=pop4(mask);
            if(v==terminal) to_t_w[u][to_t_n[u]++]=w;
            else { to_nt_v[u][to_nt_n[u]]=v; to_nt_w[u][to_nt_n[u]++]=w; }
        }
    }
    int order[70]; for(int i=0;i<S;i++) order[i]=i;
    /* sort by sum */
    for(int i=0;i<S;i++) for(int j=i+1;j<S;j++){
        int si=states[order[i]].a[0]+states[order[i]].a[1]+states[order[i]].a[2]+states[order[i]].a[3];
        int sj=states[order[j]].a[0]+states[order[j]].a[1]+states[order[j]].a[2]+states[order[j]].a[3];
        if(sj<si){ int tmp=order[i]; order[i]=order[j]; order[j]=tmp; }
    }
    int max_n=100000; int mmax=max_n+1;
    i64 *end=calloc((size_t)mmax+1,sizeof(i64));
    i64 *layers[5]; for(int i=0;i<5;i++) layers[i]=calloc((size_t)S,sizeof(i64));
    int start=idx[1][0][0][0]; layers[0][start]=1;
    for(int m=0;m<=mmax;m++){
        i64 *cur=layers[0];
        for(int oi=0;oi<S;oi++){
            int u=order[oi]; i64 val=cur[u]; if(!val) continue;
            for(int i=0;i<to_t_n[u];i++){ int nm=m+to_t_w[u][i]; if(nm<=mmax) end[nm]=(end[nm]+val)%MOD; }
            for(int i=0;i<to_nt_n[u];i++){
                int v=to_nt_v[u][i], w=to_nt_w[u][i], nm=m+w; if(nm>mmax) continue;
                if(w==0) cur[v]=(cur[v]+val)%MOD;
                else layers[w][v]=(layers[w][v]+val)%MOD;
            }
        }
        free(layers[0]);
        for(int i=0;i<4;i++) layers[i]=layers[i+1];
        layers[4]=calloc((size_t)S,sizeof(i64));
    }
    printf("%09lld\n", end[max_n+1]%MOD);
}
