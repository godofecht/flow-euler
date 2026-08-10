#include <stdint.h>
#include <stdlib.h>
typedef long long i64;
enum { MODS = 50515093, S0 = 290797, N = 10000000 };
long long pe_solve(void){
    int *lens=malloc((size_t)N*sizeof(int));
    i64 *sums=malloc((size_t)N*sizeof(i64));
    int size=0;
    i64 s=S0, prefix=0, initial_prefix_total=0;
    for(int i=0;i<N;i++){
        i64 v=s;
        prefix+=v; initial_prefix_total+=prefix;
        lens[size]=1; sums[size]=v; size++;
        while(size>=2){
            int last=size-1;
            i64 l_left=lens[last-1], t_left=sums[last-1];
            i64 l_right=lens[last], t_right=sums[last];
            i64 last_left=(t_left+l_left-1)/l_left; /* ceil */
            i64 first_right=t_right/l_right; /* floor */
            if(last_left<=first_right) break;
            /* merge */
            lens[last-1]= (int)(l_left+l_right);
            sums[last-1]=t_left+t_right;
            size--;
        }
        s=s*s%MODS;
    }
    /* final prefix sum */
    i64 final_prefix_total=0, running=0, pos_sum=0;
    for(int i=0;i<size;i++){
        i64 L=lens[i], T=sums[i];
        i64 base=T/L, rem=T%L;
        /* L-rem times base, then rem times base+1 */
        for(i64 j=0;j<L-rem;j++){ running+=base; final_prefix_total+=running; }
        for(i64 j=0;j<rem;j++){ running+=base+1; final_prefix_total+=running; }
    }
    free(lens); free(sums);
    return initial_prefix_total - final_prefix_total;
}
