lden=0;num1=1.0;num2=1.0;pprod=1.0;
R(t,K)={num=1.0;lden=0.0;p=2;for(k=1,K,num1*=(1-p^-t);num2*=(1-1/p);lden+=log(p);pprod*=p/(p-1);p=nextprime(p+1));return(num1/(num2*log(lden)))};

