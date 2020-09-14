idft(x,y,N)=
{
   local (n,m);
   for(n=1,N,y[n]=sum(m=1,N,x[m]*exp(2*Pi*I*(n-1)*(m-1)/N)));
   return(y);
}

dft(x,y,N)=
{
   local (n,m);
   for(n=1,N,y[n]=sum(m=1,N,x[m]*exp(-2*Pi*I*(n-1)*(m-1)/N)));
   return(y);
}

dif1(x,y,N)=
{
   local (n);
   for(n=1,N/2,y[n]=x[n]+x[n+N/2]);
   return(y);
}

dif2(x,y,N)=
{
   local(n);
   for(n=1,N/2,y[n]=(x[n]-x[n+N/2])*exp(2*Pi*I*(n-1)/N));
   return(y);
}

merg(xr,xi,z,N)=
{
   local(n);
   for(n=1,N,z[n]=xr[n]+I*xi[n]);
   return(z);
}

hermidft(x,y,N)=
{
   local(n);

   for(n=1,N/2,y[n]=x[n]+x[n+N/2]+exp(2*Pi*I*(n-1+N/4)/N)*(x[n]-x[n+N/2]));
   x=idft(y,x,N/2);
   return(x);
}


x=[2/16,3/16,5/16,7/16,11/16,13/16,17/16,19/16,23/16,29/16,31/16,37/16,41/16,43/16,47/16,53/16]
y=dft(x,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],16)
/*
print("y=",y)
x=idft(y,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],16)
print("Expected result=",x)
xr=dif1(y,[0,0,0,0,0,0,0,0],16)
print("xr=",xr)
xi=dif2(y,[0,0,0,0,0,0,0,0],16)
print("xi=",xi)
z=merg(xr,xi,[0,0,0,0,0,0,0,0],8);
print("z=",z)
w=idft(z,[0,0,0,0,0,0,0,0],8)
print("w=",w)
*/
z=hermidft(y,[0,0,0,0,0,0,0,0],16)