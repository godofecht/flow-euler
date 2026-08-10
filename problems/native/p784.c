#include <stdint.h>
#include <stdlib.h>
typedef long long i64;
static int bitlen_lsb(i64 x){ /* (x&-x).bit_length() */ int c=0; i64 v=x&-x; while(v){v>>=1;c++;} return c; }
long long pe_solve(void){
    i64 N=2000000;
    i64 *spf=calloc((size_t)(N+2),sizeof(i64));
    i64 *primes=malloc((size_t)(N+2)*sizeof(i64));
    i64 pc=0;
    spf[1]=1;
    for(i64 i=2;i<=N+1;i++){
        if(!spf[i]){ spf[i]=i; primes[pc++]=i; }
        for(i64 j=0;j<pc;j++){
            i64 p=primes[j];
            i64 ip=i*p;
            if(ip>N+1 || p>spf[i]) break;
            spf[ip]=p;
        }
    }
    i64 total=0;
    for(i64 r=2;r<N;r++){
        i64 kmax=N-r;
        if(kmax>r-1) kmax=r-1;
        if(kmax<=0) continue;
        i64 n_val=r*r-1;
        i64 base=2*r;
        if(kmax==1){ total+=base+1+n_val; continue; }
        i64 a=r-1,b=r+1;
        i64 fac_p[64], fac_e[64]; int nf=0;
        if(r&1){
            int ea=bitlen_lsb(a)-1, eb=bitlen_lsb(b)-1;
            int e2=ea+eb;
            a>>=ea; b>>=eb;
            if(2<=kmax){ fac_p[nf]=2; fac_e[nf]=e2; nf++; }
        }
        i64 x=a;
        while(x>1){
            i64 p=spf[x];
            if(p>kmax) break;
            int e=0; while(x>1 && spf[x]==p){ x/=p; e++; }
            fac_p[nf]=p; fac_e[nf]=e; nf++;
        }
        x=b;
        while(x>1){
            i64 p=spf[x];
            if(p>kmax) break;
            int e=0; while(x>1 && spf[x]==p){ x/=p; e++; }
            fac_p[nf]=p; fac_e[nf]=e; nf++;
        }
        i64 divs[4096]; int nd=1; divs[0]=1;
        for(int fi=0; fi<nf; fi++){
            i64 p=fac_p[fi], e=fac_e[fi];
            i64 prev_n=nd; i64 pow_p=1;
            i64 newd[4096]; int nn=0;
            for(int ee=0; ee<=e; ee++){
                for(int di=0; di<prev_n; di++){
                    i64 v=divs[di]*pow_p;
                    if(v<=kmax) newd[nn++]=v;
                }
                pow_p*=p;
                if(pow_p>kmax) break;
            }
            for(int i=0;i<nn;i++) divs[i]=newd[i];
            nd=nn;
        }
        for(int i=0;i<nd;i++) total+=base+divs[i]+n_val/divs[i];
    }
    free(spf); free(primes);
    return total;
}
