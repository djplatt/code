/* Pintz and Rusza with g=2 */
L=1;
k=19;
/* lam=3870428769./2^32; */
lam=0.3;

sig(j)=1/(k+1)*sum(m=0,k,exp(lam*cos(2*Pi*m/(k+1)))*cos(2*Pi*m*j/(k+1)));

rho(j)=-1/(k+1)*sum(m=0,k,lam*sin(2*Pi*m/(k+1))*sin(2*Pi*m*j/(k+1))*exp(lam*cos(2*Pi*m/(k+1))));

b(j)=if(j<0,0,if(j>k,0,if(j==0,sig(0),2*((k+1-j)*sig(j)-rho(j))/(k+1))));

/* this looks messy because Pari's matrices are 1 based, not 0 like normal people */


U=matrix(k,k,n,l,1/2*(b(2*(l-1)-(n-1))+b(2*(l-1)+n-1)+b(n-1-2*(l-1))));

U[1,1]=b(0);

M=U^L;

x=sum(m=1,k,M[1,m]);

print(L,"*phi(",lam,")=",log(x));

/* Heath-Brown and Puchta */

/* don't set h too big, this takes 2^h steps */

F(xi,h)=1/2^h*sum(r=0,2^h-1,exp(xi*sum(n=0,h-1,cos(2^(n-h+1)*Pi*r))));

print("log(F(",lam,",",L,"))=",log(F(lam,L)));

