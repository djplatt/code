for(n=10,2^30,fail=(1==0);forprime(p=2,sqrt(n),if(!issquarefree(n-p*p),fail=(1==1);break));if(!fail,printf("%d passed\n",n)));
quit;