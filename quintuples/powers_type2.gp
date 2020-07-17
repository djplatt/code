blim=1300000000;
count=0;

dopow(n)={
   local(lim,M,b,Mn,Mm,A,B);
   M=matrix(2,2);
   /*lim=floor((exp(log(2*blim)/n)+2)/2);*/
   lim=floor(exp(log(blim/n)/n));
   if(lim<2,printf("Nothing to be done for n=%d\n",n);return(0););
   printf("Running with r=%d for n=%d\n",lim,n);
   for(r=2,lim,
      M[1,1]=r;
      M[2,2]=r;
      ab=r*r-1;
      fordiv(ab,a,
         b=ab/a;
         M[1,2]=a;
         M[2,1]=b;
         Mn=M^n;
         A=Mn[1,2];
         B=Mn[2,1];
         if((B>A)&&(B<=blim)&&(B<A+A),
            for(m=1,n-1,
               Mm=M^m;
               count++;
               x=(Mm[1,1]-Mm[1,2])^2-1;
               if(x>=A,
                  if((x%A)==0,
                     printf("Problem with r=%d n=%d m=%d\n",r,n,m);
                     print(M," ",Mm," ",Mn)));
               if(b>a+2,
                  x=(Mm[1,1]+Mm[1,2])^2-1;
                  if(x>=A,
                     if((x%A)==0,
                        printf("Problem with r=%d n=%d m=%d\n",r,n,m);
                        print(M," ",Mm," ",Mn);return(0))))))));
   return(1);}

p=2;while(dopow(p),p=nextprime(p+1));

print("We checked a total of ",count," cases.")



