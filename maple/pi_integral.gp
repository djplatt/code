/*

Compute int(t-0.5,1.0,exp(lam^2*t^2/2)*x^t/t)

*/
\p100

X=10^23;
lam=6224003264759175.0/2^83;
lam2=lam*lam/2;
logx=log(X);

myint(a,b)=eint1(-a*logx)-eint1(-b*logx);

myintlow(a,b,m)=exp(lam2*a*a)*m;

myinthi(a,b,m)=exp(lam2*b*b)*m;

step=1/(8192*32);
resl=0;
resh=0;
forstep(st=0.5,1.0-step,step,y=st+step;m=myint(st,y);resl=resl+myintlow(st,y,m);resh=resh+myinthi(st,y,m));
print("int in[",resl,",",resh,"]");