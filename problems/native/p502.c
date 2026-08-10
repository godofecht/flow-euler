
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef int64_t i64;
#define MOD 1000000007LL
#define INV2 500000004LL

static i64 mod_pow(i64 a, i64 e, i64 mod) {
    i64 r = 1 % mod; a %= mod;
    while (e > 0) {
        if (e & 1) r = (i64)((__int128)r * a % mod);
        a = (i64)((__int128)a * a % mod);
        e >>= 1;
    }
    return r;
}
static i64 mod_inv(i64 a, i64 mod) { return mod_pow(a, mod - 2, mod); }

static void ntt(i64 *a, int n, int invert, i64 mod, i64 root) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { i64 t=a[i]; a[i]=a[j]; a[j]=t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        i64 wlen = mod_pow(root, (mod - 1) / len, mod);
        if (invert) wlen = mod_inv(wlen, mod);
        for (int i = 0; i < n; i += len) {
            i64 w = 1;
            for (int j = 0; j < len/2; j++) {
                i64 u = a[i+j];
                i64 v = (i64)((__int128)a[i+j+len/2] * w % mod);
                i64 x = u+v; if (x>=mod) x-=mod;
                i64 y = u-v; if (y<0) y+=mod;
                a[i+j]=x; a[i+j+len/2]=y;
                w = (i64)((__int128)w * wlen % mod);
            }
        }
    }
    if (invert) {
        i64 ninv = mod_inv(n, mod);
        for (int i=0;i<n;i++) a[i]=(i64)((__int128)a[i]*ninv%mod);
    }
}

static const i64 P1=998244353,G1=3,P2=1004535809,G2=3,P3=469762049,G3=3;

static i64 crt3(i64 a1,i64 a2,i64 a3){
    static i64 inv12=-1, inv123=-1; static __int128 P12;
    if(inv12<0){
        inv12=mod_inv(P1%P2,P2);
        P12=(__int128)P1*P2;
        inv123=mod_inv((i64)(P12%P3),P3);
    }
    i64 t1=(i64)((__int128)(a2-a1)%P2); if(t1<0)t1+=P2;
    t1=(i64)((__int128)t1*inv12%P2);
    __int128 x2=a1+(__int128)P1*t1;
    i64 t2=(i64)((a3-(i64)(x2%P3))%P3); if(t2<0)t2+=P3;
    t2=(i64)((__int128)t2*inv123%P3);
    __int128 x=x2+P12*t2;
    i64 r=(i64)(x%MOD); if(r<0)r+=MOD; return r;
}

static i64 *convolve_one(const i64 *a,int n1,const i64 *b,int n2,i64 mod,i64 root,int *out_n){
    int need=n1+n2-1,n=1; while(n<need)n<<=1;
    i64 *fa=calloc(n,sizeof(i64)), *fb=calloc(n,sizeof(i64));
    for(int i=0;i<n1;i++) fa[i]=a[i]%mod;
    for(int i=0;i<n2;i++) fb[i]=b[i]%mod;
    ntt(fa,n,0,mod,root); ntt(fb,n,0,mod,root);
    for(int i=0;i<n;i++) fa[i]=(i64)((__int128)fa[i]*fb[i]%mod);
    ntt(fa,n,1,mod,root);
    *out_n=need;
    i64 *res=malloc(need*sizeof(i64));
    memcpy(res,fa,need*sizeof(i64));
    free(fa); free(fb); return res;
}
static i64 *convolve_mod(const i64 *a,int n1,const i64 *b,int n2,int *out_n){
    int o1,o2,o3;
    i64 *c1=convolve_one(a,n1,b,n2,P1,G1,&o1);
    i64 *c2=convolve_one(a,n1,b,n2,P2,G2,&o2);
    i64 *c3=convolve_one(a,n1,b,n2,P3,G3,&o3);
    int need=n1+n2-1;
    i64 *res=malloc(need*sizeof(i64));
    for(int i=0;i<need;i++) res[i]=crt3(c1[i],c2[i],c3[i]);
    free(c1);free(c2);free(c3); *out_n=need; return res;
}
static i64 *poly_mul_trunc(const i64 *a,int na,const i64 *b,int nb,int n){
    int on; i64 *c=convolve_mod(a,na,b,nb,&on);
    i64 *res=calloc(n,sizeof(i64));
    for(int i=0;i<n && i<on;i++) res[i]=c[i];
    free(c); return res;
}
static i64 *inv_series(const i64 *q,int n){
    i64 *g=calloc(n,sizeof(i64));
    g[0]=mod_inv(q[0],MOD);
    int m=1;
    while(m<n){
        int m2=m*2; if(m2>n)m2=n;
        i64 *fg=poly_mul_trunc(q,m2,g,m,m2);
        for(int i=0;i<m2;i++) fg[i]=(MOD-fg[i])%MOD;
        fg[0]=(fg[0]+2)%MOD;
        i64 *ng=poly_mul_trunc(g,m,fg,m2,m2);
        free(g); free(fg); g=ng; m=m2;
    }
    return g;
}

static void build_PQ_mod(int h, int y, int deg_limit, i64 **P_out, int *Pn, i64 **Q_out, int *Qn){
    if(h<=0){
        *P_out=calloc(1,sizeof(i64)); *Pn=1;
        *Q_out=malloc(sizeof(i64)); (*Q_out)[0]=1; *Qn=1;
        return;
    }
    i64 y_mod = (y==1)?1:(MOD-1);
    int cap = deg_limit+1;
    i64 *P=calloc(cap,sizeof(i64));
    i64 *Q=calloc(cap,sizeof(i64));
    P[1]=y_mod; Q[0]=1; Q[1]=MOD-1;
    int plen= (deg_limit>=1)?2:1, qlen=2;
    if(h==1){ *P_out=P;*Pn=plen;*Q_out=Q;*Qn=qlen; return; }
    for(int step=2; step<=h; step++){
        int old_len=plen;
        int new_len=old_len+1; if(new_len>cap) new_len=cap;
        i64 *nP=calloc(new_len,sizeof(i64));
        i64 *nQ=calloc(new_len,sizeof(i64));
        nQ[0]=1;
        if(y==1){
            for(int i=1;i<new_len;i++){
                i64 s=((i-1<plen?P[i-1]:0)+(i-1<qlen?Q[i-1]:0))%MOD;
                i64 p_i=(i<old_len)?P[i]:0;
                nP[i]=(p_i+s)%MOD;
                i64 q_i=(i<old_len)?Q[i]:0;
                i64 qv=q_i-s; if(qv<0)qv+=MOD; nQ[i]=qv;
            }
        } else {
            for(int i=1;i<new_len;i++){
                i64 s=((i-1<plen?P[i-1]:0)+(i-1<qlen?Q[i-1]:0))%MOD;
                i64 p_i=(i<old_len)?P[i]:0;
                i64 v=(p_i+s)%MOD;
                nP[i]=v?MOD-v:0;
                i64 q_i=(i<old_len)?Q[i]:0;
                i64 qv=q_i-s; if(qv<0)qv+=MOD; nQ[i]=qv;
            }
        }
        free(P); free(Q); P=nP; Q=nQ; plen=qlen=new_len;
    }
    *P_out=P;*Pn=plen;*Q_out=Q;*Qn=qlen;
}

static i64 series_coeff_mod(i64 *P,int Pn,i64 *Q,int Qn,int w){
    int n=w+1;
    if(w<2048){
        i64 *f=calloc(n,sizeof(i64));
        for(int i=0;i<n;i++){
            __int128 val=(i<Pn)?P[i]:0;
            for(int j=1;j<=i && j<Qn;j++) val -= (__int128)Q[j]*f[i-j];
            f[i]=(i64)(val%MOD); if(f[i]<0)f[i]+=MOD;
        }
        i64 ans=f[w]; free(f); return ans;
    }
    i64 *Qtrunc=calloc(n,sizeof(i64));
    for(int i=0;i<n && i<Qn;i++) Qtrunc[i]=Q[i];
    i64 *invQ=inv_series(Qtrunc,n);
    i64 *Ptrunc=calloc(n,sizeof(i64));
    for(int i=0;i<n && i<Pn;i++) Ptrunc[i]=P[i];
    i64 *prod=poly_mul_trunc(Ptrunc,n,invQ,n,n);
    i64 ans=prod[w];
    free(Qtrunc); free(invQ); free(Ptrunc); free(prod);
    return ans;
}

static i64 kitamasa(i64 *init, i64 *coef, int d, i64 k){
    if(k<d) return init[k];
    i64 *pol=calloc(d,sizeof(i64)); pol[0]=1;
    i64 *base=calloc(d,sizeof(i64)); if(d>1) base[1]=1;
    i64 e=k;
    while(e){
        if(e&1){
            /* mul_reduce pol*base */
            i64 *tmp=calloc(2*d-1,sizeof(i64));
            for(int i=0;i<d;i++) if(pol[i])
                for(int j=0;j<d;j++) tmp[i+j]=(tmp[i+j]+(i64)((__int128)pol[i]*base[j]%MOD))%MOD;
            for(int i=2*d-2;i>=d;i--){
                i64 t=tmp[i];
                if(t) for(int j=0;j<d;j++) tmp[i-1-j]=(tmp[i-1-j]+(i64)((__int128)t*coef[j]%MOD))%MOD;
            }
            memcpy(pol,tmp,d*sizeof(i64)); free(tmp);
        }
        {
            i64 *tmp=calloc(2*d-1,sizeof(i64));
            for(int i=0;i<d;i++) if(base[i])
                for(int j=0;j<d;j++) tmp[i+j]=(tmp[i+j]+(i64)((__int128)base[i]*base[j]%MOD))%MOD;
            for(int i=2*d-2;i>=d;i--){
                i64 t=tmp[i];
                if(t) for(int j=0;j<d;j++) tmp[i-1-j]=(tmp[i-1-j]+(i64)((__int128)t*coef[j]%MOD))%MOD;
            }
            memcpy(base,tmp,d*sizeof(i64)); free(tmp);
        }
        e>>=1;
    }
    i64 ans=0;
    for(int i=0;i<d;i++) ans=(ans+(i64)((__int128)pol[i]*init[i]%MOD))%MOD;
    free(pol); free(base); return ans;
}

static i64 coeff_large_w_small_h(i64 w,int h,int y){
    i64 *P,*Q; int Pn,Qn;
    build_PQ_mod(h,y,h,&P,&Pn,&Q,&Qn);
    i64 *f=malloc((h+1)*sizeof(i64));
    for(int n=0;n<=h;n++){
        __int128 val=(n<Pn)?P[n]:0;
        for(int j=1;j<=n && j<Qn;j++) val -= (__int128)Q[j]*f[n-j];
        f[n]=(i64)(val%MOD); if(f[n]<0)f[n]+=MOD;
    }
    if(w<=h){ i64 ans=f[w]; free(P);free(Q);free(f); return ans; }
    i64 *coef=malloc(h*sizeof(i64));
    for(int i=1;i<=h;i++) coef[i-1]=(MOD-(i<Qn?Q[i]:0))%MOD;
    i64 *init=malloc(h*sizeof(i64));
    for(int i=0;i<h;i++) init[i]=f[i+1];
    i64 ans=kitamasa(init,coef,h,w-1);
    free(P);free(Q);free(f);free(coef);free(init);
    return ans;
}

static i64 *poly_mul_small(const i64 *a,int na,const i64 *b,int nb,int limit,int *out_n){
    int need=na+nb-1; if(need>limit) need=limit;
    i64 *res=calloc(need,sizeof(i64));
    for(int i=0;i<na;i++) if(a[i]){
        for(int j=0;j<nb;j++){
            int ij=i+j; if(ij>=limit) break;
            res[ij]=(res[ij]+(i64)((__int128)a[i]*b[j]%MOD))%MOD;
        }
    }
    *out_n=need; return res;
}
static i64 *poly_add_small(const i64 *a,int na,const i64 *b,int nb,int limit){
    i64 *res=calloc(limit,sizeof(i64));
    for(int i=0;i<limit;i++){
        i64 av=i<na?a[i]:0, bv=i<nb?b[i]:0;
        res[i]=(av+bv)%MOD;
    }
    return res;
}

typedef struct { i64 *p; int n; } Poly;
typedef struct { Poly m[2][2]; } Mat;

static Poly Pmake(i64 *p,int n){ Poly r={p,n}; return r; }
static void Pfree(Poly p){ free(p.p); }

static Mat mat_mul(Mat A, Mat B, int limit){
    Mat C;
    for(int i=0;i<2;i++) for(int j=0;j<2;j++){
        int n1,n2;
        i64 *t1=poly_mul_small(A.m[i][0].p,A.m[i][0].n,B.m[0][j].p,B.m[0][j].n,limit,&n1);
        i64 *t2=poly_mul_small(A.m[i][1].p,A.m[i][1].n,B.m[1][j].p,B.m[1][j].n,limit,&n2);
        i64 *s=poly_add_small(t1,n1,t2,n2,limit);
        free(t1); free(t2);
        C.m[i][j]=Pmake(s,limit);
    }
    return C;
}
static void mat_free(Mat M){
    for(int i=0;i<2;i++) for(int j=0;j<2;j++) Pfree(M.m[i][j]);
}
static Mat mat_pow(Mat M, i64 e, int limit){
    Mat I;
    for(int i=0;i<2;i++) for(int j=0;j<2;j++){
        i64 *p=calloc(limit,sizeof(i64));
        if(i==j) p[0]=1;
        I.m[i][j]=Pmake(p,limit);
    }
    Mat base=M;
    int base_owned=0;
    while(e){
        if(e&1){
            Mat nI=mat_mul(I,base,limit);
            mat_free(I); I=nI;
        }
        Mat nB=mat_mul(base,base,limit);
        if(base_owned) mat_free(base);
        base=nB; base_owned=1;
        e>>=1;
    }
    if(base_owned) mat_free(base);
    return I;
}

static i64 coeff_large_h_small_w(int w, i64 h, int y){
    if(h<=0) return 0;
    int limit=w+1;
    i64 y_mod=(y==1)?1:(MOD-1);
    i64 *P=calloc(limit,sizeof(i64));
    i64 *Q=calloc(limit,sizeof(i64));
    if(w>=1){ P[1]=y_mod; Q[1]=MOD-1; }
    Q[0]=1;
    if(h==1){
        i64 ans=series_coeff_mod(P,limit,Q,limit,w);
        free(P); free(Q); return ans;
    }
    /* Build transfer matrix A */
    Mat A;
    {
        i64 *y_one_plus_x=calloc(limit,sizeof(i64));
        y_one_plus_x[0]=y_mod;
        if(limit>=2) y_one_plus_x[1]=y_mod;
        i64 *y_x=calloc(limit,sizeof(i64));
        if(limit>=2) y_x[1]=y_mod;
        i64 *minus_x=calloc(limit,sizeof(i64));
        if(limit>=2) minus_x[1]=MOD-1;
        i64 *one_minus_x=calloc(limit,sizeof(i64));
        one_minus_x[0]=1;
        if(limit>=2) one_minus_x[1]=MOD-1;
        A.m[0][0]=Pmake(y_one_plus_x,limit);
        A.m[0][1]=Pmake(y_x,limit);
        A.m[1][0]=Pmake(minus_x,limit);
        A.m[1][1]=Pmake(one_minus_x,limit);
    }
    Mat Mp=mat_pow(A,h-1,limit);
    mat_free(A);
    /* apply to (P,Q) */
    int n1,n2,n3,n4;
    i64 *t1=poly_mul_small(Mp.m[0][0].p,limit,P,limit,limit,&n1);
    i64 *t2=poly_mul_small(Mp.m[0][1].p,limit,Q,limit,limit,&n2);
    i64 *Pn=poly_add_small(t1,n1,t2,n2,limit);
    free(t1); free(t2);
    t1=poly_mul_small(Mp.m[1][0].p,limit,P,limit,limit,&n3);
    t2=poly_mul_small(Mp.m[1][1].p,limit,Q,limit,limit,&n4);
    i64 *Qn=poly_add_small(t1,n3,t2,n4,limit);
    free(t1); free(t2);
    mat_free(Mp);
    free(P); free(Q);
    i64 ans=series_coeff_mod(Pn,limit,Qn,limit,w);
    free(Pn); free(Qn);
    return ans;
}

static i64 coeff_C_mod(i64 w, i64 h, int y){
    if(h<=0) return 0;
    if(w<=200 && h>2000) return coeff_large_h_small_w((int)w,h,y);
    if(h<=200 && w>2000) return coeff_large_w_small_h(w,(int)h,y);
    i64 *P,*Q; int Pn,Qn;
    int deg=(int)w; if(deg<0)deg=0;
    build_PQ_mod((int)h,y,deg,&P,&Pn,&Q,&Qn);
    i64 ans=series_coeff_mod(P,Pn,Q,Qn,(int)w);
    free(P); free(Q);
    return ans;
}
static i64 E_leq_mod(i64 w,i64 h){
    if(h<=0) return 0;
    i64 T=coeff_C_mod(w,h,1);
    i64 S=coeff_C_mod(w,h,-1);
    return (i64)((__int128)(T+S)%MOD * INV2 % MOD);
}
static i64 F_mod(i64 w,i64 h){
    i64 a=E_leq_mod(w,h);
    i64 b=E_leq_mod(w,h-1);
    i64 r=a-b; if(r<0)r+=MOD; return r;
}

long long pe502_answer(void){
    i64 a=F_mod(1000000000000LL,100);
    i64 b=F_mod(10000,10000);
    i64 c=F_mod(100,1000000000000LL);
    return (a+b+c)%MOD;
}
