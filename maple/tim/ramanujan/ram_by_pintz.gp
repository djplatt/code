H=1.2e12;
RHbound=solve(lx=10,100,4.92*sqrt(exp(lx)/lx)-H);
A=195.6;B=2.014;C=1.963; /* valid from exp(4000) */

e0(x,lx)={local(R,X);if(lx<RHbound,exp(-lx/2)*lx^2/(8*Pi),if(lx<4000,X=sqrt(lx/5.69693);sqrt(8*X/Pi)*exp(-X),R=5.573412;A*(log(x)/R)^B*exp(-C*sqrt(log(x)/R))))};


a(x)={local(lx);lx=log(x);lx^5*e0(x,lx)};


eMa(lx)=72+2*M+(2*M+132)/lx+(4*M+288)/lx^2+(12*M+576)/lx^3+48*M/lx^4+M^2/lx^5;
ema(lx)=206+m+364/lx+381/lx^2+238/lx^3+97/lx^4+30/lx^5+8/lx^6;

th(t)={local(res);res=0;forprime(p=2,t,res+=log(p));res};
C2=intnum(t=2,149,(th(t)-t+t*e0(t))/(t*log(t)^2));
C3=2*sum(k=1,5,k!/log(2)^(k+1));

/* do a search */
for(t=4040260,4040270,{lxa=t/1000.0;
xa=exp(lxa);
axa=a(xa);
C1=intnum(t=log(3),lxa,exp(t)*a(exp(t))/t^7)+(2-log(2))*log(3)^5; /* last bit 2->3 */
M=120+axa+C1*lxa^6/xa+(720+axa)*(1/lxa+7*2^8/lxa^2+7*lxa^6/(sqrt(xa)*log(2)^8));
m=120-axa-(C1+C2+C3)*lxa^6/xa+(720-axa)*(1/lxa+7*2^8/lxa^2+7*lxa^6/(sqrt(xa)*log(2)^8));
lxdash=eMa(lxa)-ema(lxa);
print("log x' = ",lxdash, " log xa+1 = ",lxa+1," so we can take x = exp(",max(lxa+1,lxdash),")");});


/* this is a good value */
lxa=4040.266;
xa=exp(lxa);
axa=a(xa);
C1=intnum(t=log(3),lxa,exp(t)*a(exp(t))/t^7)+(2-log(2))*log(3)^5;
M=120+axa+C1*lxa^6/xa+(720+axa)*(1/lxa+7*2^8/lxa^2+7*lxa^6/(sqrt(xa)*log(2)^8));
m=120-axa-(C1+C3)*lxa^6/xa+(720-axa)*(1/lxa+7*2^8/lxa^2+7*lxa^6/(sqrt(xa)*log(2)^8));
lxdash=eMa(lxa)-ema(lxa);
print("log x' = ",lxdash, " log xa+1 = ",lxa+1," so we can take x = exp(",max(lxa+1,lxdash),")");



