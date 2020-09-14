
/*
The original e0 function to use when log(x) is less than 35
*/
e00(x)={local(R,X);R=6.455;X=sqrt(log(x)/R);sqrt(8*X/(17*Pi))*exp(-X)}

/*
A revised one using Buethe/FK bounds
*/
e0(x)={local(R,X,lx);lx=log(x);if(lx>=100,1.75185e-10,if(lx>=90,1.79330e-10,if(lx>=80,1.84848e-10,if(lx>=70,1.9191e-10,if(lx>=65,1.96865e-10,if(lx>=55,2.88434e-10,if(lx>=50,1.16465e-9,if(lx>=45,1.11742e-8,if(lx>=40,8.6347e-8,if(lx>=35,7.4457e-7,e00(x)))))))))))};
/*
e0(x)={local(R,X,lx);lx=log(x);if(lx>=100,1.75185e-10,if(lx>=90,1.79330e-10,if(lx>=80,1.84848e-10,if(lx>=70,1.9191e-10,if(lx>=65,1.96865e-10,if(lx>=55,2.88434e-10,if(lx>=50,1.16465e-9,if(lx>=45,4.82e-9,if(lx>=44,6.62e-9,if(lx>=43,9.58e-9,if(lx>=42,1.45e-8,if(lx>=41,2.26e-8,if(lx>=40,3.63e-8,if(lx>=39,5.26e-8,if(lx>=35,3.08e-7,e00(x))))))))))))))))};
*/
a(x)=log(x)^5*e0(x);

pip(x)=(x+x*e0(x))/log(x)+intnum(t=log(2),log(x),exp(t)/t^2)+intnum(t=log(1.39e17),log(x),exp(t)*e0(exp(t))/t^2);

th(t)={local(res);res=0;forprime(p=2,t,res+=log(p));res};

X0=149;
th_X0=intnum(t=2,X0,th(t)/(t*log(t)^2));

pim(x)=(x-x*e0(x))/log(x)+th_X0+intnum(t=log(X0),log(x),exp(t)/t^2)-intnum(t=log(X0),log(x),exp(t)*e0(exp(t))/t^2);

g(x)=pip(x)^2-exp(1)*x/log(x)*pim(x/exp(1));


eMa(lx)=72+2*M+(2*M+132)/lx+(4*M+288)/lx^2+(12*M+576)/lx^3+48*M/lx^4+M^2/lx^5;
ema(lx)=206+m+364/lx+381/lx^2+238/lx^3+97/lx^4+30/lx^5+8/lx^6;

C2=intnum(t=2,149,(th(t)-t+t*e0(t))/(t*log(t)^2));
C3=2*sum(k=1,5,k!/log(2)^(k+1));

/*
Test log(x) in ranges of width 0.1
Bound the troublesome integral (C4) by taking it all the way
to x_a+0.1
*/
for(t=350,489,{lxa=t/10.;
C4=intnum(t=log(2),lxa+.1,exp(t)/t^7);
xa=exp(lxa);
axa=a(xa);
C1=intnum(t=log(2),lxa,exp(t)*a(exp(t))/t^7);
M=120+axa+C1*lxa^6/xa+(720+axa)*(C4*lxa^6/xa);
m=120-axa-(C1+C2+C3)*lxa^6/xa+(720-axa)*(C4*lxa^6/xa);
lxdash=eMa(lxa)-ema(lxa);
printf("log x' = %5.4f log xa+1 = %5.4f so we can take x = exp(%5.4f)\n",lxdash,lxa+1,max(lxa+1,lxdash));});
/*
Same again but with bigger steps to log(x)=132. I reckon we can stop sooner
but whatever.
*/
for(t=49,132,{lxa=t;
C4=intnum(t=log(2),lxa+1,exp(t)/t^7);
xa=exp(lxa);
axa=a(xa);
C1=intnum(t=log(2),lxa,exp(t)*a(exp(t))/t^7);
M=120+axa+C1*lxa^6/xa+(720+axa)*(C4*lxa^6/xa);
m=120-axa-(C1+C2+C3)*lxa^6/xa+(720-axa)*(C4*lxa^6/xa);
lxdash=eMa(lxa)-ema(lxa);
printf("log x' = %5.4f log xa+1 = %5.4f so we can take x = exp(%5.4f)\n",lxdash,lxa+1,max(lxa+1,lxdash));});

