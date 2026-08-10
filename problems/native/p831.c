#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
enum { DEG = 5 };
static void poly_mul(mpz_t *res, mpz_t *a, mpz_t *b) {
    mpz_t tmp[DEG+1];
    for(int i=0;i<=DEG;i++){ mpz_init(tmp[i]); mpz_set_ui(tmp[i],0); }
    for(int i=0;i<=DEG;i++){
        if(mpz_sgn(a[i])==0) continue;
        for(int j=0;j<=DEG-i;j++){
            if(mpz_sgn(b[j])==0) continue;
            mpz_addmul(tmp[i+j], a[i], b[j]);
        }
    }
    for(int i=0;i<=DEG;i++){ mpz_set(res[i], tmp[i]); mpz_clear(tmp[i]); }
}
static void poly_pow(mpz_t *res, mpz_t *base, int exp) {
    mpz_t b[DEG+1], r[DEG+1];
    for(int i=0;i<=DEG;i++){ mpz_init_set(b[i], base[i]); mpz_init(r[i]); mpz_set_ui(r[i], i==0?1:0); }
    while(exp>0){
        if(exp&1) poly_mul(r,r,b);
        exp>>=1;
        if(exp) poly_mul(b,b,b);
    }
    for(int i=0;i<=DEG;i++){ mpz_set(res[i], r[i]); mpz_clear(b[i]); mpz_clear(r[i]); }
}
void pe_solve_print(void){
    int m=142857;
    mpz_t Q[DEG+1], Qm[DEG+1], A[DEG+1], AQm[DEG+1];
    int Qc[6]={1,3,5,5,3,1};
    int Ac[6]={1,5,10,10,5,1}; /* C(5,k) */
    for(int i=0;i<=DEG;i++){
        mpz_init_set_ui(Q[i], Qc[i]);
        mpz_init(Qm[i]);
        mpz_init_set_ui(A[i], Ac[i]);
        mpz_init(AQm[i]);
    }
    poly_pow(Qm, Q, m);
    poly_mul(AQm, A, Qm);
    /* convert AQm[5] to base 7, first 10 digits */
    mpz_t c, rem; mpz_init_set(c, AQm[DEG]); mpz_init(rem);
    if(mpz_sgn(c)==0){ printf("0000000000\n"); return; }
    char digits[256]; int len=0;
    while(mpz_sgn(c)>0){
        mpz_fdiv_qr_ui(c, rem, c, 7);
        digits[len++] = (char)('0'+mpz_get_ui(rem));
    }
    /* digits are least-significant first; reverse for MSB */
    char out[32]; int o=0;
    for(int i=len-1;i>=0 && o<10;i--) out[o++]=digits[i];
    while(o<10) out[o++]='0';
    out[10]=0;
    printf("%s\n", out);
    for(int i=0;i<=DEG;i++){ mpz_clear(Q[i]); mpz_clear(Qm[i]); mpz_clear(A[i]); mpz_clear(AQm[i]); }
    mpz_clear(c); mpz_clear(rem);
}
