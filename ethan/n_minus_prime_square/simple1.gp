for(n=2^30+1,2^30+2^28,fail=(1==0);forprime(p=2,sqrt(n),if(!issquarefree(n-p*p),fail=(1==1);break));if(!fail,printf("%d passed\n",n)));
quit;