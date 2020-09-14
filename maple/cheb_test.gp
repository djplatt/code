s0=0;s1=0;s2=0;s3=0;x0=10^6;p=2;while(p<2*x0,s0++;s1+=p-x0;s2+=(p-x0)^2;s3+=(p-x0)^3;print(p," ",s1," ",s2," ",s3);p=nextprime(p+1));
quit;