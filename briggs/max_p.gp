/*
Given an eps compute the exponents of the primes of the respective 
CA number
*/

/* set working precision */
\p100

/* ap as defined in Theorem 2.1 (i) */
logp(x,p)=log(x)/log(p);
ap1(eps,p)=logp((p^(1+eps)-1)/(p^eps-1),p);
ap(eps,p)=floor(ap1(eps,p))-1;

/*
gen_ca(eps)={
   local(last_p,last_ap,p,this_ap,np);
   if(eps<1e-102,printf("solve in gen_ca won't work with eps below about 1e-102.\n");return(0));
   last_p=2;last_ap=ap(eps,2);p=3;this_ap=ap(eps,3);
   if(last_ap==0,printf("1\n");return(0));
   while(this_ap!=last_ap,printf("%d^%d ",last_p,last_ap);last_ap=this_ap;last_p=p;p=nextprime(p+1);this_ap=ap(eps,p);if(this_ap<=1,break));
   while(this_ap>0,np=precprime(solve(pp=2,1e100,ap1(eps,pp)-this_ap-1));if(last_p!=np,printf("[%d..%d]^%d ",last_p,np,this_ap),printf("%d^%d ",p,this_ap));this_ap--;last_p=nextprime(np+1));
printf("\n")}
*/

gen_ca(eps)={
   local(last_p,last_ap,p,this_ap,np);
   if(eps<1e-102,printf("solve in gen_ca won't work with eps below about 1e-102.\n");return(0));
   last_p=2;last_ap=ap(eps,2);p=3;this_ap=ap(eps,3);
   if(last_ap==0,printf("1\n");return(0));
   while(this_ap!=last_ap,printf("%d %d %d\n",last_p,last_p,last_ap);last_ap=this_ap;last_p=p;p=nextprime(p+1);this_ap=ap(eps,p);if(this_ap<=1,break));
   while(this_ap>0,np=precprime(solve(pp=2,1e100,ap1(eps,pp)-this_ap-1));if(last_p!=np,printf("%d %d %d\n",last_p,np,this_ap),printf("%d %d %d\n",p,p,this_ap));this_ap--;last_p=nextprime(np+1));
printf("\n")}

gen_ca(1e-12);
quit;



