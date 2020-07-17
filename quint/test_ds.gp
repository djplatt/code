allocatemem(8000000000);

X=readvec("solutions.lst");

print("We have ",length(X)/4," potential solutions.");

dp(a,b,c)={r=sqrtint(a*b+1);s=sqrtint(a*c+1);t=sqrtint(b*c+1);return(a+b+c+2*(r*s*t+a*b*c));}

border(a,b,d)=0.0033*sqrt(a*d)/(4*b)<log(4.001*a*b^2*d)*log(1.299*sqrt(a*b)*d/(b-a))/(log(4*b*d)*log(0.1053*d/(b*(b-a)^2)));

/*p=1;for(i=1,length(X)/4,a=X[p];p++;b=X[p];p++;c=X[p];p++;d=X[p];p++;if(d!=dp(a,b,c),print(a," ",b," ",c)));*/

fails=0;p=1;for(i=1,length(X)/4,a=X[p];p++;b=X[p];p++;c=X[p];p++;d=X[p];p++;if((b>=1.45*a)&&!border(a,b,d),fails++;));print("We had ",fails," failures.");

