allocatemem(2^30)

\p100

N=2^18

h=176.5

t0=1000000.0

one_over_A=41.0/1024.0

g(t,t0,h,k)=(-2*Pi*t*I)^k*exp(lngamma(0.25+(t+t0)/2*I)+Pi*(t+t0)/4-t^2/(2*h^2))

g_k=-1

make_g_vec(t0,h,k)=
{
   g_vec=vector(N,n,g((n-N/2-1)*one_over_A,t0,h,k));
   g_k=k;
}

G(n,k)=
{
   local(om,res,m,om_m);

   if(g_k!=k,make_g_vec(t0,h,k));

   om=exp(-2.0*Pi*I*n/N);
   print(om);
   om_m=1.0;
   res=0;
   for(m=0,N-1,
      res=res+g_vec[m+1]*om_m;
      om_m=om_m*om;);
   res=res*one_over_A/k!;
   if((n%2)==1,return(-res),return(res));
}