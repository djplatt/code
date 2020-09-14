N_t(q,t)=
{
   return(t/2/Pi*log(q*t/2/Pi/exp(1)));
}

num_prims1(p,n)=
{
   local (phi,res);

   if(n==1,return(p-2));
   phi=p^(n-1)*(p-1);
   res=phi-phi/p;
   return(res);
}

num_prims(N)=
{
   local(facs,res);

   facs=factorint(N);
   res=1;
   for(i=1,matsize(facs)[1],res=res*num_prims1(facs[i,1],facs[i,2]));
   return(res);
}

num_real_prims(N)=
{
   local(facs);

   facs=factorint(N);
   if(facs[1,1]==2,if(facs[1,2]==1,return(0),if(facs[1,2]==2,res=1,if(facs[1,2]==3,res=2,return(0)))),if(facs[1,2]==1,res=1,return(0)));
   for(i=2,matsize(facs)[1],if(facs[i,2]>1,return(0)));
   return(res);
}

num_pairs(N)=
{
   local(np,nr);

   np=num_prims(N);
   nr=num_real_prims(N);
   return((np-nr)/2+nr);
}

num_zeros(q,t)=
{
   return(t/(2*Pi)*log(q*t/(2*Pi*exp(1))));
}