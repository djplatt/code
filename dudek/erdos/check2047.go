for(n=1,2047,if(n%4!=1,p=2;p2=p*p;while(p2<n,if(issquarefree(n-p2),print("Success ",n," ",p);break,p=nextprime(p+1);p2=p*p));if(p2>=n,print("Failure at n=",n))));

quit;
