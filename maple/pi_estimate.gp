\p200

lim1=10000;
lim=10400000;

p=2;lres=0;hres=0;psum=0;lp=log(p);

while(p<lim1,np=nextprime(p+1);lnp=log(np);psum+=lp;w=psum*intnum(y=p,np,1/(y*log(y)*log(y)));lres+=w;p=np;lp=lnp)

print(lres);

hres=lres;

while(p<lim,np=nextprime(p+1);lnp=log(np);psum+=lp;lres+=psum*(np-p)/(np*lnp*lnp);hres+=psum*(np-p)/(p*lp*lp);lp=lnp;p=np);
/*
*intnum(y=p,np,1/(y*(log(y))^2));p=np)
*/
print("int theta(y)/y/log^2 y from 2 to ",lim,"=[",lres,",",hres,"]");

