R=5.573412;
myroundup(x,n)={local(e);e=floor(log(x)/log(10))-n;return(ceil(x/10^(e))*10^(e))};
myrounddown(x,n)={local(e);e=floor(log(x)/log(10))-n;return(floor(x/10^(e))*10^(e))};

delta(lx)=2/sqrt(R*lx);
T(lx)=exp(2*sqrt(lx/R));
C(lx)=2*sqrt(1-8*delta(lx)/3);
B(lx)=2+delta(lx);
A1(lx)=exp(2)/2*2.2*2^(4+2*delta(lx))/R^(2+delta(lx));
A2(lx)=exp(2)*6*2^3*R^(-1.5)*lx^(1.5-2-delta(lx))*exp(-2*sqrt(lx/R))*exp(2*sqrt((1-8*delta(lx)/3)*lx/R));
A3(lx)=2/Pi*(lx/R)^(-1-delta(lx))*exp(-2*sqrt(lx/R))*exp(2*sqrt((1-8*delta(lx)/3)*lx/R));

A(lx)=A1(lx)+A2(lx)+A3(lx);

ABC(lx)=printf("%d & %4.3f & %4.3f & %4.3f \\\\\n",lx,myroundup(A(lx),3),myroundup(B(lx),3),myroundup(C(lx),3));

/*for(n=0,18,ABC(1000+500*n));*/

/* Version 9 */
AA(lx)=12*(4)^(2+delta(lx));
BB(lx)=2+delta(lx);
CC(lx)=2*sqrt(1-8*delta(lx)/3)/sqrt(1+1/(lx-1));

ABC1(lx)=printf("$%d$ & %4.1f & %4.3f & %4.3f \\\\\n",lx,myroundup(AA(lx),4),myroundup(BB(lx),3),myrounddown(CC(lx),3));

for(n=0,9,ABC1(1000+1000*n));
