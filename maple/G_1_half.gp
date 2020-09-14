x=10^22-273747*2^32+1;

logx=log(x);

lam=7699934548453755/2^79;
l2=lam*lam/2;

resl=0;
resh=0;

int(N,A,B)=
{
   local (step,C1,C2,ei);
   resl=0;
   resh=0;
   step=(B-A)/N;
   C1=A;
   C2=A+step;
   while(C2<=B,ei=eint1(-C2*logx)-eint1(-C1*logx);resl=resl-exp(l2*C1*C1)*ei;resh=resh-exp(l2*C2*C2)*ei;C1=C2;C2=C2+step);
   print("int in [",resl,",",resh,"]");
}

