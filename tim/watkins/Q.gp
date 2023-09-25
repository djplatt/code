 d=-1832763;

alpha(a,d)=log(a)+log(8*Pi*exp(-Euler))-1/2*log(d);

H(a,b,d)=
{
   k=sqrt(d)/2/a;
   return((2*cos(Pi*b/a)/sqrt(k)/exp(2*Pi*k)+2*abs(cos(Pi*b/a))/25/sqrt(k)/exp(2*Pi*k)));
}

do_sums(a,b,c,d)=
{
  printf("%d %d %d %d %f %f\n",a,b,c,d,alpha(a,-d)/sqrt(a),H(a,b,-d)/sqrt(a));
}

check_me(a,b,c,d)={
   if(((-a<b)&&(b<=a)&&(a<c))||((0<=b)&&(b<=a)&&(a==c)),do_sums(a,b,c,d));
}

for(b=0,sqrt(-d/3),ac4=b*b-d;if(ac4%4==0,ac=ac4/4;fordiv(ac,a,c=ac/a;check_me(a,b,c,d);check_me(a,-b,c,d))));

quit;



