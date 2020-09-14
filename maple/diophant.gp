
forprime(p=2,100,print(p," ",1.0*sum(x=0,p-1,sum(y=0,p-1,sum(z=0,p-1,if((x^3+y^3-z^3+33)%p==0,1))))/p^2));
/*

N=7*13*11;
H=10000;
W=33;
for(x=0,N-1,for(y=0,N-1,for(z=0,N-1,if((x^3+y^3-z^3+33)%N==0,print("./timamod ",W," ",H," ",N," ",x," ",y," ",z)))));
quit;
*/