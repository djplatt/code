gb(n)={c=0;if(isprime(n-4),c++);forprime(p=ceil(n/3),n-6,forprime(q=ceil((n-p)/2),p,if(isprime(n-p-q),c++)));print(n," ",c);}

/* for(n=3,5000,gb(2*n+1));

quit;
*/

S(x,alpha)={tot=0;forprime(p=2,floor(x),tot+=exp(2*Pi*I*alpha*p));tot;}

N=50000.0;for(a=0,floor(N/2),print(a/N," ",abs(intnum(al=a/N,(a+1)/N,S(97,al)^3*exp(-2*Pi*I*97*al)))));

/*for(n=0,1000,al=n/2000.0;print(al," ",log(abs(S(97,al)^3*exp(-2*Pi*I*97*al)))));*/
quit;


