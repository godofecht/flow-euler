#include <stdint.h>
#include <string.h>
typedef long long i64;
enum { MAXL = 25, MAXS = 9*MAXL };
static int is_prime_arr[MAXS+1];
static i64 counts[MAXL+1][MAXS+1];
static void init(void){
    memset(is_prime_arr,1,sizeof is_prime_arr);
    is_prime_arr[0]=is_prime_arr[1]=0;
    for(int p=2;p*p<=MAXS;p++) if(is_prime_arr[p]) for(int j=p*p;j<=MAXS;j+=p) is_prime_arr[j]=0;
    memset(counts,0,sizeof counts);
    counts[0][0]=1;
    for(int L=1;L<=MAXL;L++){
        for(int s=0;s<=9*L;s++){
            i64 total=0;
            for(int d=0;d<=9 && d<=s;d++) total+=counts[L-1][s-d];
            counts[L][s]=total;
        }
    }
}
static i64 count_len(int length){
    if(length<=0) return 0;
    if(length==1){ i64 t=0; for(int d=1;d<=9;d++) if(is_prime_arr[d]) t++; return t; }
    i64 total=0; int rem=length-1;
    for(int first=1;first<=9;first++)
        for(int ts=0;ts<=9*rem;ts++) if(is_prime_arr[first+ts]) total+=counts[rem][ts];
    return total;
}
static i64 kth_len(int length, i64 k){
    int digits[MAXL]; int prefix_sum=0;
    for(int pos=0;pos<length;pos++){
        int rem=length-pos-1; int start=pos==0?1:0;
        for(int d=start;d<=9;d++){
            i64 cnt=0; int base=prefix_sum+d;
            for(int ts=0;ts<=9*rem;ts++) if(is_prime_arr[base+ts]) cnt+=counts[rem][ts];
            if(k>cnt) k-=cnt;
            else { digits[pos]=d; prefix_sum+=d; break; }
        }
    }
    i64 ans=0; for(int i=0;i<length;i++) ans=ans*10+digits[i];
    return ans;
}
long long pe_solve(void){
    init();
    i64 remaining=10000000000000000LL;
    for(int length=1;length<=MAXL;length++){
        i64 cnt=count_len(length);
        if(remaining>cnt) remaining-=cnt;
        else return kth_len(length, remaining);
    }
    return -1;
}
