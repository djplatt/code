\p 30
myroundup(x,n)=floor(x/10^floor(log(x)/log(10)-n+1)+1)*10.^floor(log(x)/log(10)-n+1);
myrounddown(x,n)=ceil(x/10^floor(log(x)/log(10)-n+1)-1)*10.^floor(log(x)/log(10)-n+1);
B=myrounddown(solve(x=1e10,1e30,4.92*sqrt(x/log(x))-3e12),4);
print("B=",B);


pn=29996208012611;

print("pn=",pn);

R=5.573412;
BB=1.52;
C=myrounddown(1.89/sqrt(R),3);
A=myroundup(411.5/R^BB,3);
printf("A=%f B=%f C=%f\n",A,BB,C);

lX0=2000;
e55=myroundup(1.388e-10+1.4262/sqrt(exp(55)),4);
print("eps(55)=",e55);
alp=myroundup(e55*(1+log(B))/log(B),4);
print("alpha=",alp);
C1a=myroundup(alp*(log(lX0)-log(log(B))),4);
print("C1a=",C1a);
C1b=myroundup(A*2/C*exp(-C*sqrt(lX0))*(6/C^3+6*sqrt(lX0)/C^2+3*lX0/C+lX0^1.5),4);
print("C1b=",C1b);
C1=myroundup(C1a+C1b,4);
print("C1=",C1);

gB(t,pn,B)={local(lpn,spn);lpn=log(pn);spn=sqrt(pn);exp(2/pn)*lpn*exp(1.02/(pn-1)/lpn+lpn/spn/8/Pi+C1+((lpn+3)*sqrt(B)-(log(B)+3)*sqrt(pn))/(4*Pi*sqrt(pn*B)))/(zeta(t)*log(pn-spn*lpn^2/(8*Pi)))};

ginf(t,pn)={local(lpn);lpn=log(pn);log(pn)*exp(2/pn+1.02/(pn-1)/lpn+1/6/lpn^3+5/8/lpn^4)/(zeta(t)*log(pn-1.4262*sqrt(pn)-1.338e-10*pn))};

print("gB(20)=",gB(20,pn,B));

print("ginf(20)=",ginf(20,B));

print("gB(21)=",gB(21,pn,B));

print("ginf(21)=",ginf(21,B));