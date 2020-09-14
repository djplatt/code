/* 
   |theta(x)-x|/x <= e0(x,log(x)) for x>=2 
   trivial bound to 599
   RH bound above
*/

e0(x,lx)={local(R,X);if(x<599,(2-log(2))/2,exp(-lx/2)*lx^2/(8*Pi))};

a(lx,x)=lx^5*e0(x,lx);

eMa(lx)={local(M);M=MM(lx);return(72+2*M+(2*M+132)/lx+(4*M+288)/lx^2+(12*M+576)/lx^3+48*M/lx^4+M^2/lx^5)};
ema(lx)={local(m);m=mm(lx);return(206+m+364/lx+381/lx^2+238/lx^3+97/lx^4+30/lx^5+8/lx^6)};

MM(lxa)={local(axa,xa);xa=exp(lxa);axa=a(lxa,xa);return(120+axa+C1(lxa)*lxa^6/xa+(720+axa)*(1/lxa+7*2^8/lxa^2+7*lxa^6/(sqrt(xa)*log(2)^8)))};
mm(lxa)={local(axa,xa);xa=exp(lxa);axa=a(lxa,xa);return(120-axa-(C1(lxa)+C2)*lxa^6/xa+(-axa)*(1/lxa+7*2^8/lxa^2+7*lxa^6/(sqrt(xa)*log(2)^8)))};

C1(lxa)=intnum(t=log(2),lxa,exp(t)*a(t,exp(t))/t^7);

/* not needed */
th(t)={local(res);res=0;forprime(p=2,t,res+=log(p));res};

C2=2*sum(k=1,5,k!/log(2)^(k+1));

/* do a search */

for(t=1360,1380,\
lxa=t/10.0;\
lxdash=eMa(lxa)-ema(lxa);\
print("log x' = ",lxdash, " log xa+1 = ",lxa+1," so we can take x = exp(",max(lxa+1,lxdash),")"));



