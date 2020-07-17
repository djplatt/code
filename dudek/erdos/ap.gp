/*

Compute some better constants for Adrian's attack on
Erdos conjecture n=p^2 + square free, n>9 n % 4 <neq 1

*/
 
num_prims1(p,n)=
{
   local (phi,res);

   phi=p^(n-1)*(p-1);
   res=phi-ceil(phi/p);
   return(res);
}

/* how many primitive characters does N have */
num_prims(N)=
{
   local(facs,res);

   facs=factorint(N);
   res=1;
   for(i=1,matsize(facs)[1],res=res*num_prims1(facs[i,1],facs[i,2]));
   return(res);
}

two_pi=2*Pi;

C1=9.14; /* R+R 4.2. Theorem 3.6.3 would be better */
C2=0.9185;
C3=5.512; /* McCurley from Ramare+Rumely lemma 4.1.1 */
R=9.645908801;


/* How far did I test GRH to for modulus k
   In fact we only used H=1000
*/

round_h(h0)=ceil(h0/10)*10;

H(k)=if(k>100000,print("modulus too large, need GRH to 1000. Exiting");quit,if(k==1,30610046000,if(k%2==1,round_h(max(10^8/k,37500000/k+200)),round_h(max(10^8/k,75000000/k+200)))));

/* R+R 4.3.4 */

A(m,delta)=((1+(1+delta)^(m+1))/delta)^m;
/*
This is the definition in R+R 4.1.4, but we use Lemma 4.2.1
*/
/*
At(h,k,m,x0)=intnum(t=h,[+1],(log(k*t/two_pi)/Pi+C2/t)*exp(-log(x0)/(R*log(k*t/C1)))/t^(m+1));
*/

/* 4.2.1 */
hm(t,k,m)=1/(Pi*(m-1)*t^(m-1))*(log(k*t/two_pi)+1/(m-1))+C2/(m*t^m);

/* Lemma 4.2.1 */

At(h,k,m,x0)=hm(h,k,m)/h*exp(-log(x0)/(R*log(k*h/C1)));

/* defined in Lemma 4.1.3 */

Bt(h,k,m,x0)=1/h^(m+1)*exp(-log(x0)/(R*log(k*h/C1)))*2*(C2*log(k*h)+C3);

Ct(h,k,m)=1/(Pi*m*h^m)*(log(k*h/two_pi)+1/m);

Dt(h,k,m)=1/h^(m+1)*(2*C2*log(k*h)+2*C3+C2/(m+1));

/* Et(chi) = sum 1/|rho| |Im rho|<=H
   but we use Lemma 4.1.2 as a bound.
*/ 
Et(h,k)=(1/two_pi*log(h)^2+log(k/two_pi)/Pi*log(h)+C2+2*(log(k/(two_pi*exp(1)))/Pi+C2*log(k)+C3));

/* 4.3.1 */
f(k)={local(res);res=0.0;fordiv(k,d,if(isprime(d),res+=1/(d-1)));return(res)};

/* 4.3.2 */
Rt(k,x0)=eulerphi(k)/x0*((f(k)+0.5)*log(x0)+4*log(k)+13.4);

/* contribution from a primitive character of modulus k Th 5.1.1 */
per_ABCD_chi(h,k,m,x0)=At(h,k,m,x0)+Bt(h,k,m,x0)+(Ct(h,k,m)+Dt(h,k,m))/sqrt(x0);

/* contribution from all characters modulus k. sum over all the primitive characters modulo d where d|k
*/
all_ABCD_chi(h,k,m,x0)={local(res);res=0.0;fordiv(k,d,if(d%4!=2,res+=per_ABCD_chi(h,d,m,x0)*num_prims(d)));return(res)};

/* summing Et over all characters modulo k */
all_E_chi(h,k)={local(res);res=0.0;fordiv(k,d,if(d%4!=2,res+=Et(h,d)*num_prims(d)));return(res)};

/* find an epsilon for psi Theorem 5.1.1 */
eps(h,k,m,x0,delta)={local(t1,t2,t3,t4);t1=1/2*A(m,delta)*all_ABCD_chi(h,k,m,x0);t2=(1+m*delta/2)/sqrt(x0)*all_E_chi(h,k);t3=m*delta/2;t4=Rt(k,x0);/*print(t1," ",t2," ",t3," ",t4," ",t1+t2+t3+t4);*/return(t1+t2+t3+t4)};

/* find an epsilon for theta just below table 1, page 420*/
ept(h,k,m,x0,delta)=eps(h,k,m,x0,delta)+eulerphi(k)*(1.0012/x0^0.5+3/x0^(2/3));
\p200

/* it turns out max for x<=10^10, |theta(x;q^2,l)-x/phi(q^2)|/sqrt(x) happens with x=7, l=7 for all q in [17,19,...97]
The following eq works for x in sqrt(2.5e14) to 10^10
and says that |theta-x/phi| < eq x/phi
*/
myeq(q)=(log(7)-7/(q*(q-1)))/sqrt(7)*(q*(q-1))/(2.5e14)^0.25;

myround(x)=ceil(x*100000)/100000.;

/* the primes q in R+R tables already covered */
ps=[3,5,7,11,13];

/* data for q^2 from R+R with table 2 figure for 25 corrected */

table1=[0.003228,0.012214,0.017015,0.031939,0.042497];
table2=[1.108042,0.821891,0.744132,0.711433,0.718525];

/* sum the contribution from these initial primes */
printf("Computing eq for moduli q^2 with q in [3..13] to apply for x>2.5e14^0.5\n");
res=0.0;for(n=1,5,q=ps[n];eq=table2[n]*q*(q-1)/(2.5e14)^0.25;eq=max(eq,table1[n]);eq=myround(eq);printf("eq for q=%d is %7.5f\n",q,eq);res+=2*(1+eq)/(q*(q-1)));

printf("Initial error from primes in [3...13]=%7.5f\n",myround(res));

forprime(q=17,97,eq=myround(myeq(q));q2=q*q;h=1000;delta=2*exp(1)/h;best=ept(h,q2,2,10^10,delta);bestm=1;for(m=3,30,t=ept(h,q2,m,10^10,delta);if(t<best,best=t;bestm=m));best=myround(best);printf("Best for q= %d at m= %d was ",q,bestm);rbest=max(best,eq);printf("%7.5f vs %7.5f using eq=%7.5f\n",best,eq,rbest);res+=2*(1+rbest)/(q*(q-1)));printf("New error for primes in [3...97]=%7.5f\n",myround(res));

quit;

