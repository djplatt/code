X=readvec("$HOME/data/dudek/fails.lst");

for(n=1,length(X),p=2;while(!issquarefree(X[n]-p*p),p=nextprime(p+1));print(X[n]," ",p," ",X[n]%8));


quit;
