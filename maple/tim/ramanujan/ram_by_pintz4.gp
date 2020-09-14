/* height to which we assume RH */
H=2.5e12;
/* corresponding x to which we can use RH bound on theta */ 
RHbound=floor(solve(lx=10,100,4.92*sqrt(exp(lx)/lx)-H)); /* Buthe */
print("Taking RH to be true to exp(",RHbound,")");
A=356.4;B=1.52;C=1.89; /* valid from exp(4000) - DJP/TST Pintz*/

/* 
   |theta(x)-x|/x <= e0(x,log(x)) for x>=2 
   trivial bound to 599
   RH bound to Buthe bound on partial RH
   Trudgian Th 1 bound to exp(4000)
   DJP/TST Pintz bound beyond
*/

e0(x,lx)={local(R,X,A,B,C);if(x<599,(2-log(2))/2,if(lx<RHbound,exp(-lx/2)*lx^2/(8*Pi),if(lx<1169,X=sqrt(lx/6.455);sqrt(8*X/(17*Pi))*exp(-X),R=5.573412;B=1.52;C=1.89;if(lx<2000,A=462.0,if(lx<3000,A=411.5,if(lx<4000,A=379.7,A=356.4)));A*exp(-C*sqrt(lx/R))*(lx/R)^B)))};

a(lx,x)=lx^5*e0(x,lx);

eMa(lx)={local(M);M=MM(lx-1);return(72+2*M+(2*M+132)/lx+(4*M+288)/lx^2+(12*M+576)/lx^3+48*M/lx^4+M^2/lx^5)};
ema(lx)={local(m);m=mm(lx-1);return(206+m+364/lx+381/lx^2+238/lx^3+97/lx^4+30/lx^5+8/lx^6)};

MM(lxa)={local(axa,xa);xa=exp(lxa);axa=a(lxa,xa);return(120+a(lxa+1,xa*exp(1))+C1(lxa)*lxa^6/xa+(720+axa)*(1/lxa+7*2^8/lxa^2+7*lxa^6/(sqrt(xa)*log(2)^8)))};

mm(lxa)={local(axa,xa);xa=exp(lxa);axa=a(lxa,xa);return(120-a(lxa+1,xa*exp(1))-(C2(lxa)+C3)*lxa^6/xa-axa*(1/lxa+7*2^8/lxa^2+7*lxa^6/(sqrt(xa)*log(2)^8)))};

C1(lxa)=intnum(t=log(2),lxa,(720+exp(t)*a(t,exp(t)))/t^7);
C2(lxa)=intnum(t=log(2),lxa,(720-exp(t)*a(t,exp(t)))/t^7);

/* not needed */
th(t)={local(res);res=0;forprime(p=2,t,res+=log(p));res};

C3=2*sum(k=1,5,k!/log(2)^(k+1));
print("C3=",C3);

/* do a search */
/*
for(t=4003150,4003160,\
lxa=t/1000.0;\
lxdash=eMa(lxa)-ema(lxa);\
print("log x' = ",lxdash, " log xa+1 = ",lxa+1," so we can take x = exp(",max(lxa+1,lxdash),")"));
*/

print("Ramaujan is true above exp(",solve(lx=2000,5000,eMa(lx)-ema(lx)-lx),")");


