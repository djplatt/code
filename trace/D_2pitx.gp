hhat(t)=if(t<a,2*(33*a^5-30*a^3*t^2+15*a*t^4-5*t^5),if(t<2*a,51*a^5+75*a^4*t-210*a^3*t^2+150*a^2*t^3-45*a*t^4+5*t^5,(3*a-t)^5))/(120*a^6);

L(t)=
{
	local(d,l,Lval,ff);
	d=quaddisc(t^2-4);
	l=round(sqrt((t^2-4)/d));
	Lval=2/sqrt(d)*quadregulator(d)*qfbclassno(d);
	ff=factor(l);
	Lval *= prod(i=1,matsize(ff)[1],
		1+(ff[i,1]-kronecker(d,ff[i,1]))*(ff[i,1]^ff[i,2]-1)/(ff[i,1]-1)
	);
	Lval/l
}

H(maxt)=
{
   local(res,tterm,t,lv);
   t=3;
   res=0;
   tterm=log((t+sqrt(t^2-4))/2)/Pi;
   while(tterm<maxt,
      res+=L(t)*hhat(tterm);
      t++;
      tterm=log((t+sqrt(t^2-4))/2)/Pi;);
   print("Last t used was ",t-1);
   return(res/Pi);
};

H1_pow(i,maxt)=
{
   local(res,p,lp,lterm);
   p=2;
   n=p^i;
   lp=log(p);
   lterm=i*lp/Pi;
   while(lterm<maxt,
      res+=lp/n*hhat(lterm);
      p=nextprime(p+1);
      lp=log(p);
      lterm=i*lp/Pi);
   print("last p^",i," used was ",precprime(p-1));
   return(res);
}

H1(maxt)=
{
   local(res,i,h1);
   i=1;
   res=0;
   h1=H1_pow(1,maxt);
   while(h1>0,res+=h1;i++;h1=H1_pow(i,maxt));
   return(res/(Pi));
}

a=1/4;

print("a=",a);
print("D=",H(3*a)+H1(3*a));

quit;
