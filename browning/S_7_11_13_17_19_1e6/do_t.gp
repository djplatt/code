H=1000000;
lH=log(H);
first_p=7;
last_p=19;
N=1;
forprime(p=first_p,last_p,N*=p^(floor(lH/log(p))));
fordiv(N,n,if(n<=H,print(n)));
quit;
