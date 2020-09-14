lim=10000;p=2;th=log(2);res=0;
forprime(q=3,lim,res+=th*intnum(t=p,q,1/(t*log(t)^2));th+=log(q);p=q);
res+=th*intnum(t=p,lim,1/(t*log(t)^2));
