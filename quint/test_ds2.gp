allocatemem(8000000000);

X=readvec("$HOME/data/quint/solutions.log");

print("We have ",length(X)/4," potential solutions.");

dp(a,b,c)={r=sqrtint(a*b+1);s=sqrtint(a*c+1);t=sqrtint(b*c+1);return(a+b+c+2*(r*s*t+a*b*c));}

border(a,b,d,g3)=g3*sqrt(a*d)/(4*b)<log(4.001*a*b^2*d)*log(1.299*sqrt(a*b)*d/(b-a))/(log(4*b*d)*log(0.1053*d/(b*(b-a)^2)));

/*p=1;for(i=1,length(X)/4,a=X[p];p++;b=X[p];p++;c=X[p];p++;d=X[p];p++;if(d!=dp(a,b,c),print(a," ",b," ",c)));*/


G1(a,b,c)={local(t1,t2,t3,res,n);t1=log(2.001^2*b*c);t2=log(1.994*sqrt(a)/sqrt(b));t3=log(1.994^2*a*c);n=ceil(solve(nn=2,1000000,(nn+1)/nn-(t1+t2/nn)/t3));return((t1+t2/n)/t3);}

G2(a,b,g1)={local(r,g12);r=b/a*1.;g12=g1^2;if(r>g12,return(r),repl++;return(g12));}

G3(a,c,g1,g2)=1/g1*sqrt(1+1/(a*c))*(sqrt((g2*a*c+1)/(a*c+1))-g1);

/*
pass=0;fails=0;p=1;for(i=1,length(X)/4,a=X[p];p++;b=X[p];p++;c=X[p];p++;d=X[p];p++;g1=G1(a,b,d);g2=G2(a,b,g1);g3=G3(a,d,g1,g2);if(!border(a,b,d,g3),fails++,pass++;print(pass,": ",fails," : ( ",a, " , ",b, " , ",c, " , ",d," ) is still a possible.")));print("We had ",fails," failures.");print("We had ",pass," passes.");


quit;
*/
