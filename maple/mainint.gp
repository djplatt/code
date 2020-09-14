\p200
x=10.^23;
lam=6224003264759175.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/8.0;

lam2=lam*lam/2;
logx=log(x);

step=1/512.0;

res=0.0;

forstep(st=0.5,1.0-step,step,res=res+intnum(t=st,st+step,exp(lam2*t*t+t*logx)/t));

print("int=",res);
