#include <stdint.h>
typedef long long i64;
static int digit_sum(i64 n){i64 s=0; while(n){s+=n%10;n/=10;} return (int)s;}
static int dfs(i64 num, i64 target){
    if(target<0||target>num) return 0;
    if(num==target) return 1;
    if(num==0) return target==0;
    if((num-target)%9!=0) return 0;
    if(target<120 && digit_sum(num)>target) return 0;
    for(i64 pow10=10; pow10<=num; pow10*=10){
        i64 suffix=num%pow10;
        if(suffix>target) break;
        i64 prefix=num/pow10;
        if(prefix==0) break;
        if(dfs(prefix, target-suffix)) return 1;
    }
    return 0;
}
static int is_s_root(i64 root){
    i64 sq=root*root;
    if(sq<10) return 0;
    for(i64 pow10=10; pow10<=sq; pow10*=10){
        i64 suffix=sq%pow10;
        if(suffix>root) break;
        i64 prefix=sq/pow10;
        if(prefix==0) break;
        if(dfs(prefix, root-suffix)) return 1;
    }
    return 0;
}
long long pe_solve(void){
    i64 limit=1000000000000LL;
    i64 max_root=1000000; /* sqrt(1e12) */
    i64 total=0;
    for(i64 base=0; base<=1; base++){
        i64 start = base==0 ? 9 : 1;
        for(i64 r=start; r<=max_root; r+=9){
            if(r<4) continue;
            i64 sq=r*r;
            if(sq<=limit && is_s_root(r)) total+=sq;
        }
    }
    return total;
}
