/*
e0(x)={local(R,X);R=6.455;X=sqrt(log(x)/R);sqrt(8*X/(17*Pi))*exp(-X)}
*/
e0(x)={local(R,X);R=5.573412;X=sqrt(log(x)/R);sqrt(8*X/Pi)*exp(-X)}

a(x)=log(x)^5*e0(x);

pip(x)=(x+x*e0(x))/log(x)+intnum(t=log(2),log(x),exp(t)/t^2)+intnum(t=log(1.39e17),log(x),e0(exp(t))/t^2);

th(t)={local(res);res=0;forprime(p=2,t,res+=log(p));res};

X0=149;
th_X0=intnum(t=2,X0,th(t)/(t*log(t)^2));

pim(x)=(x-x*e0(x))/log(x)+th_X0+intnum(t=log(X0),log(x),exp(t)/t^2)-intnum(t=log(X0),log(x),e0(exp(t))/t^2);

g(x)=pip(x)^2-exp(1)*x/log(x)*pim(x/exp(1));


eMa(lx)=72+2*M+(2*M+132)/lx+(4*M+288)/lx^2+(12*M+576)/lx^3+48*M/lx^4+M^2/lx^5;
ema(lx)=206+m+364/lx+381/lx^2+238/lx^3+97/lx^4+30/lx^5+8/lx^6;

C2=intnum(t=2,149,(th(t)-t+t*e0(t))/(t*log(t)^2));
C3=2*sum(k=1,5,k!/log(2)^(k+1));

/*
for(t=8800030,8800050,{lxa=t/1000.0;
xa=exp(lxa);
axa=a(xa);
C1=intnum(t=log(2),lxa,exp(t)*a(exp(t))/t^7);
M=120+axa+C1*lxa^6/xa+(720+axa)*(1/lxa+7*2^8/lxa^2+7*lxa^6/(sqrt(xa)*log(2)^8));
m=120-axa-(C1+C2+C3)*lxa^6/xa+(720-axa)*(1/lxa+7*2^8/lxa^2+7*lxa^6/(sqrt(xa)*log(2)^8));
lxdash=eMa(lxa)-ema(lxa);
print("log x' = ",lxdash, " log xa+1 = ",lxa+1," so we can take x = exp(",max(lxa+1,lxdash),")");});
*/

lxa=8800.037;
xa=exp(lxa);
axa=a(xa);
C1=intnum(t=log(2),lxa,exp(t)*a(exp(t))/t^7);

M=120+axa+C1*lxa^6/xa+(720+axa)*(1/lxa+7*2^8/lxa^2+7*lxa^6/(sqrt(xa)*log(2)^8));
m=120-axa-(C1+C2+C3)*lxa^6/xa+(720-axa)*(1/lxa+7*2^8/lxa^2+7*lxa^6/(sqrt(xa)*log(2)^8));
lxdash=eMa(lxa)-ema(lxa);
print("log x' = ",lxdash, " log xa+1 = ",lxa+1," so we can take x = exp(",max(lxa+1,lxdash),")");



