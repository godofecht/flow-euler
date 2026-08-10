#include <stdint.h>
long long pe_solve(void){
    const int n=1000000, k=10;
    double m=2.0*k;
    double mu=0,s2=0;
    double M=(double)k*(n+1);
    double total=0,comp=0;
    double kf=(double)k;
    for(int i=0;i<n;i++){
        double Ey=mu+kf;
        double Ey2=s2+2*kf*mu+kf*kf;
        double denom=M*(M-1);
        double c1=m*(M-m)/denom;
        double c2=m*(m-1)/denom;
        double Eb2=c1*Ey+c2*Ey2;
        double y=Eb2-comp;
        double t=total+y;
        comp=(t-total)-y;
        total=t;
        double alpha=(M-m)/M;
        double mu_next=alpha*Ey;
        double coeff_y2=1.0-2.0*m/M+c2;
        double s2_next=coeff_y2*Ey2+c1*Ey;
        mu=mu_next; s2=s2_next; M-=kf;
    }
    return (long long)(total+0.5);
}
