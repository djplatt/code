/*
Lehman's theorem with new s_1

See Saouter & Demichel 2010

*/

/* A is height to which RH is known
   T is height to which we have zeros
*/


check_params(al,et,om,A,T)=
{

/*   if(om-et<=25.57,print("omega-eta must be >= 25.57");return(0)); */
   if(om-et<=43.613,print("omega-eta must be >= 43.613");return(0));
   if(4*A/om>al,print("4A/omega must be <= alpha");return(0));
   if(al>A*A,print("Alpha must be <= A^2");return(0));
   if(2*A/al>et,print("2A/alpha must be <= eta");return(0));
   if(et>om/2,print("eta must be <= omega/2");return(0));
   return(1);
}


S1(et,om)=2/(om-et)+10.04/(om-et)^2+log(2)*(om+et)*exp(-(om-et)/2)+2/log(2)*(om+et)*exp(-(om-et)/6);

S2(al,et)=2*exp(-al*et^2/2)/sqrt(2*Pi*al)/et;

S3(al,et)=0.08*sqrt(al)*exp(-al*et^2/2);

S4(al,T)=exp(-T^2/2/al)*(al/Pi/T^2*log(T/2/Pi)+8*log(T)/T+4*al/T^3);

S5(et,om)=0.05/(om-et);

S6(al,et,om,A,T)=A*log(A)*exp(-A^2/2/al+(om+et)/2)*(4/sqrt(al)+15*et);

R(al,et,om,A,T)=if(check_params(al,et,om,A,T),S1(et,om)+S2(al,et)+S3(al,et)+S4(al,T)+S5(et,om)+S6(al,et,om,A,T),0);

K(y,al)=sqrt(al/2/Pi)*exp(-al*y*y/2);








/* A new improved version using Dusart 2010 */

/*new_S1(et,om)=2/(om-et)+9.336/(om-et)^2+log(2)*(om+et)*exp(-(om-et)/2)+2/log(2)*(om+et)*exp(-(om-et)/6);*/
/* Esteki Introduced Nov 2013 */
new_S1(et,om)=2/(om-et)+8/(om-et)^2+(7.5976*8)/(om-et)^3+log(2)*(om+et)*exp(-(om-et)/2)+2/log(2)*(om+et)*exp(-(om-et)/6);
/* Esteki True for x>10^150 (see Trudgian email 26/11/13*/
new_S1(et,om)=2/(om-et)+(2.022*4)/(om-et)^2+log(2)*(om+et)*exp(-(om-et)/2)+2/log(2)*(om+et)*exp(-(om-et)/6);


/* not needed at all */

new_S2(al,et)=0;

/* Better approximation to li(exp(rho u)) */

new_S3(al,et,om)=K(et,al)*(0.2+0.016/(om-et));

/* =(om-et)^-2*2*2*sum gamma^-3 First 2 is for +ve and -ve gamma */
new_S5(et,om)=0.0032/(om-et)^2;

/* from Zegowitz */
new_S6(al,et,om,A,T)=A*log(A)*exp(-A^2/2/al+(om+et)/2)*(3.2/sqrt(al)+14.4*et);


new_R(al,et,om,A,T)=if(check_params(al,et,om,A,T),new_S1(et,om)+new_S2(al,et)+new_S3(al,et,om)+S4(al,T)+new_S5(et,om)+new_S6(al,et,om,A,T),0);

tails(al,et,et0,om)=(et-et0)*K(et0,al)*(0.51*(exp((om-et0)/2)/(om-et)^2+exp((om+et)/2)/(om+et0)^2)+1.80141*((om+et)*exp(-(om+et0)/2)+(om-et0)*exp(-(om-et)/2)));
/*
print("Checking Saouter and Demichel original error terms");
A=6.85e7;
al=6e12;
T=13046000.0
om=727.951335792;
et=2*A/al;
print("omega=",om);
print("eta=",et);
print("alpha=",al);
print("A=",A);
print("T=",T);
print("R=",R(al,et,om,A,T));
sm=-1.0028892566;
print("Sum over zeros came to ",sm);
print("Margin=",1+R(al,et,om,A,T)+sm);
*/

print("");
print("Checking Saouter and Demichel new error terms");
A=6.85e7;
al=6e12;
T=13046000.0
om=727.951335792;
et=0.00002283333334;
print("omega=",om);
print("eta=",et);
print("alpha=",al);
print("A=",A);
print("T=",T);
print("new_S1=",new_S1(et,om));
print(2/(om-et)+9.336/(om-et)^2);
print("S2=",S2(al,et));
print("New S3=",new_S3(al,et,om));
print("S4=",S4(al,T));
print("new_S5=",new_S5(et,om));
print("New S6=",new_S6(al,et,om,A,T));

print("R=",new_R(al,et,om,A,T));

sm=-1.0028892566;
print("Sum over zeros came to ",sm);
print("Margin=",1+new_R(al,et,om,A,T)+sm);
et0=et/2.077;
print("tails=",tails(al,et,et0,om));
print("With new et0=",1+new_R(al,et,om,A,T)+sm+tails(al,et,et0,om));
print("min pi-li=",-(1+new_R(al,et,om,A,T)+sm+tails(al,et,et0,om))*exp((om-et0)/2)/(om-et0));

/*
print("");
print("Checking Left of Saouter and Demichel new error terms");
A=6.85e7;
al=6e12;
T=13046000.0
om=727.951335766;
et=0.00002283333334;
print("omega=",om);
print("eta=",et);
print("alpha=",al);
print("A=",A);
print("T=",T);
print("new_S1=",new_S1(et,om));
print(2/(om-et)+9.336/(om-et)^2);
print("S2=",S2(al,et));
print("New S3=",new_S3(al,et,om));
print("S4=",S4(al,T));
print("new_S5=",new_S5(et,om));
print("New S6=",new_S6(al,et,om,A,T));

print("R=",new_R(al,et,om,A,T));

sm=-1.0027779809716;
print("Sum over zeros came to ",sm);
print("Margin=",1+new_R(al,et,om,A,T)+sm);
*/


/*
print("");
print("Checking Right of Saouter and Demichel new error terms");
A=6.85e7;
al=6e12;
T=13046000.0
om=727.951335796;
et=0.00002283333334;
print("omega=",om);
print("eta=",et);
print("alpha=",al);
print("A=",A);
print("T=",T);
print("new_S1=",new_S1(et,om));
print(2/(om-et)+9.336/(om-et)^2);
print("S2=",S2(al,et));
print("New S3=",new_S3(al,et,om));
print("S4=",S4(al,T));
print("new_S5=",new_S5(et,om));
print("New S6=",new_S6(al,et,om,A,T));

print("R=",new_R(al,et,om,A,T));

sm=-1.002906481004;
print("Sum over zeros came to ",sm);
print("Margin=",1+new_R(al,et,om,A,T)+sm);
print("With new et0=",1+new_R(al,et,om,A,T)+sm+tails(al,et,et/2.078,om));
print("om+et0=",om+et/2.078);
*/


print("");
print("Checking Improvement of Saouter and Demichel new error terms, new alpha,A,T");
f=1.1;
A=6.85e7*f;
al=6e12*f;
T=13046000+2100000*6;
om=727.951335780;
et=0.00002283333334;
print("omega=",om);
print("eta=",et);
print("alpha=",al);
print("A=",A);
print("T=",T);
print("new_S1=",new_S1(et,om));
print(2/(om-et)+9.336/(om-et)^2);
print("S2=",S2(al,et));
print("New S3=",new_S3(al,et,om));
print("S4=",S4(al,T));
print("new_S5=",new_S5(et,om));
print("New S6=",new_S6(al,et,om,A,T));

print("R=",new_R(al,et,om,A,T));

sm=-1.00282390423749793;
print("Sum over zeros came to ",sm);
print("Margin=",1+new_R(al,et,om,A,T)+sm);
f=solve(f=1,4,1+new_R(al,et,om,A,T)+sm+tails(al,et,et/f,om));
f=floor(f*100000);
print("Reducing to eta_0 by a factor of ",f,"/100000");
f=f/100000.0;
print("Tail integrals bounded by ",tails(al,et,et/f,om));
print("With new et0, Margin=",1+new_R(al,et,om,A,T)+sm+tails(al,et,et/f,om));

print("om-et0=",floor((om-et/f)*10000000000),"/10000000000");
print("om+et0=",ceil((om+et/f)*10000000000),"/10000000000");


/*
print("");
print("Checking Sharpened Saouter and Demichel new error terms");
A=3e10;
T=378446000;
om=727.95134682+727.95132478;
om=om/2;
et=727.95134682-om;
al=ceil(2*A/et);
print("omega=",om);
print("eta=",et);
print("alpha=",al);
print("A=",A);
print("T=",T);
print("R=",new_R(al,et,om,A,T));
sm=0.0;
print("Margin=",1+new_R(al,et,om,A,T)+sm);

print("new_S1=",new_S1(et,om));
print("S2=",S2(al,et));
print("New S3=",new_S3(al,et,om));
print("S4=",S4(al,T));
print("new_S5=",new_S5(et,om));
print("New S6=",new_S6(al,et,om,A,T));

/*
A=30610046000;
T=3001346000;
om=727.95132478+727.95134681;
om=om/2;
al=A^2/om;
al=al/4;
et=2.*A/al;

print("Saouter and Demichel Sharpened error terms:- R=",new_R(al,et,om,A,T));
print("new_S1=",new_S1(et,om));
print("S2=",S2(al,et));
print("New S3=",new_S3(al,et,om));
print("S4=",S4(al,T));
print("new_S5=",new_S5(et,om));
print("S6=",S6(al,et,om,A,T));
print("omega=",om);
print("eta=",et);
print("alpha=",al);
print("A=",A);
print("T=",T);

A=30610046000;
T=3001346000;
om=727.95133307663053;
al=2^57;
et=2.*A/al;

print("Demichel region, new error terms:- R=",new_R(al,et,om,A,T));
print("new_S1=",new_S1(et,om));
print("S2=",S2(al,et));
print("New S3=",new_S3(al,et,om));
print("S4=",S4(al,T));
print("new_S5=",new_S5(et,om));
print("New S6=",new_S6(al,et,om,A,T));
print("eta=",et);

A=30610046000;
T=3001346000;
om=727.95133307663053;
al=A^2/om;
al=al/4;
et=2.*A/al;
print("Demichel region, new error terms, tighter eta:- R=",new_R(al,et,om,A,T));
print("new_S1=",new_S1(et,om));
print("S2=",S2(al,et));
print("New S3=",new_S3(al,et,om));
print("S4=",S4(al,T));
print("new_S5=",new_S5(et,om));
print("New S6=",new_S6(al,et,om,A,T));
print("eta=",et);
print("alpha=",al);
print("omega=",om);
print("T=",T);
print("A=",A);
*/

roundn(x,n)={
local(count);
count=0;
while(x<0.5,x=x+x;count--);
while(x>1,x=x/2;count++);
x*=2^n;
x=round(x);
return(x*2^(count-n)*1.);
};
/*
om=727.95133598;
A=3.0610046e10;
print("Setting up for omega=",om," and A=",A);
al=solve(al=A^2/2/om,A^2/om,new_S6(al,0,om,A,A)-1e-15);
print("alpha was ",al);
al=round(al);
print("Alpha set to ",al);
et=roundn(2*A/al,20);
print("Eta set to ",et);
T=solve(T=A/100,A,S4(al,T)-1e-10);
T=ceil((T-446000)/2100000)*2100000+446000
print("T set to ",T);
print("last file is zeros_",T-21000000);

print("new_S1=",new_S1(et,om));
print("S2=",S2(al,et));
print("New S3=",new_S3(al,et,om));
print("S4=",S4(al,T));
print("new_S5=",new_S5(et,om));
print("New S6=",new_S6(al,et,om,A,T));

print("R=",new_R(al,et,om,A,T));

li(x)=-eint1(-log(x));
riemann_R(x)=sum(n=2,floor(log(x)/log(2)),-moebius(n)/n*li(x^(1/n)));

pi_upb1(x)=li(x)-riemann_R(x)/exp(1)*(log(log(log(x)))+1-exp(1));
pi_upb2(x)=x/log(x)*(1+2/log(x)+2.334/(log(x))^2);
*/

print("Version of October 2012");
print("Checking Improvement of Saouter and Demichel new error terms, new alpha,A,T,omega,eta");

om=727.95133598;
A=3.0610046e10;

al=1.1533087226142278600000000000000000000000000000e18;
T=6.989246000000000000e+09;
et=5.3082146678207209333777427673339843750000000000e-8;
print("omega=",om);
print("eta=",et);
print("alpha=",al);
print("A=",A);
print("T=",T);
print("new_S1=",new_S1(et,om));
print("new_S2=",new_S2(al,et));
print("New S3=",new_S3(al,et,om));
print("S4=",S4(al,T));
print("new_S5=",new_S5(et,om));
print("New S6=",new_S6(al,et,om,A,T));

print("R=",new_R(al,et,om,A,T));

sm=-1.004318016142537621045743;
print("Sum over zeros came to ",sm);
print("Margin=",1+new_R(al,et,om,A,T)+sm);
f=solve(f=1,4,1+new_R(al,et,om,A,T)+sm+tails(al,et,et/f,om));
f=floor(f*100000);
print("Reducing to eta_0 by a factor of ",f,"/100000");
f=f/100000.0;
print("Tail integrals bounded by ",tails(al,et,et/f,om));
m=-1-new_R(al,et,om,A,T)-sm-tails(al,et,et/f,om);
print("With new et0, Margin=",m);

print("om-et0=",floor((om-et/f)*10000000000),"/10000000000");
print("om+et0=",ceil((om+et/f)*10000000000),"/10000000000");

print(m*exp(floor((om-et/f)*10000000000)/10000000000/2)/ceil((om+et/f)*10000000000)*10000000000);

print("Version of October 2012");
print("Checking Improvement of Saouter and Demichel new error terms, new alpha,A,T,omega,eta");

om=727.95133588;
A=3.0610046e10;

al=1.1533087226142278600000000000000000000000000000e18;
T=6.989246000000000000e+09;
et=5.3082146678207209333777427673339843750000000000e-8;
print("omega=",om);
print("eta=",et);
print("alpha=",al);
print("A=",A);
print("T=",T);
print("new_S1=",new_S1(et,om));
print("new_S2=",new_S2(al,et));
print("New S3=",new_S3(al,et,om));
print("S4=",S4(al,T));
print("new_S5=",new_S5(et,om));
print("New S6=",new_S6(al,et,om,A,T));

print("R=",new_R(al,et,om,A,T));

sm=-1.00318521254458;
print("Sum over zeros came to ",sm);
print("Margin=",+1+new_R(al,et,om,A,T)+sm);
f=solve(f=1,4,1+new_R(al,et,om,A,T)+sm+tails(al,et,et/f,om));
f=floor(f*100000);
print("Reducing to eta_0 by a factor of ",f,"/100000");
f=f/100000.0;
print("Tail integrals bounded by ",tails(al,et,et/f,om));
m=-1-new_R(al,et,om,A,T)-sm-tails(al,et,et/f,om);
print("With new et0, Margin=",m);

print("om-et0=",floor((om-et/f)*10000000000),"/10000000000");
print("om+et0=",ceil((om+et/f)*10000000000),"/10000000000");

print(m*exp(floor((om-et/f)*10000000000)/10000000000/2)/ceil((om+et/f)*10000000000)*10000000000);

print("Version of November 2012");
print("Checking best to date");

om=727.951333078;
A=3.0610046e10;

al=1.1533087226142278600000000000000000000000000000e18;
T=6.989246000000000000e+09;
et=5.3082146678207209333777427673339843750000000000e-8;
print("omega=",om);
print("eta=",et);
print("alpha=",al);
print("A=",A);
print("T=",T);
print("new_S1=",new_S1(et,om));
print("new_S2=",new_S2(al,et));
print("New S3=",new_S3(al,et,om));
print("S4=",S4(al,T));
print("new_S5=",new_S5(et,om));
print("New S6=",new_S6(al,et,om,A,T));

print("R=",new_R(al,et,om,A,T));

sm=-1.00280277617;
print("Sum over zeros came to ",sm);
print("Margin=",1+new_R(al,et,om,A,T)+sm);
f=solve(f=1,4,1+new_R(al,et,om,A,T)+sm+tails(al,et,et/f,om));
f=floor(f*100000);
print("Reducing to eta_0 by a factor of ",f,"/100000");
f=f/100000.0;
print("Tail integrals bounded by ",tails(al,et,et/f,om));
m=-1-new_R(al,et,om,A,T)-sm-tails(al,et,et/f,om);
print("With new et0, Margin=",m);

print("om-et0=",floor((om-et/f)*10000000000),"/10000000000");
print("om+et0=",ceil((om+et/f)*10000000000),"/10000000000");

print(m*exp(floor((om-et/f)*10000000000)/10000000000/2)/ceil((om+et/f)*10000000000)*10000000000);

print("Version of November 2013");
print("Searching for a new record!");

om=727.95133290;
A=3.0610046e10;

al=A^2/1000;
/*1.1533087226142278600000000000000000000000000000e18;*/
T=1.1155646000e+10;
et=2.0*A/al;
/*5.3082146678207209333777427673339843750000000000e-8;*/
print("omega=",om);
print("eta=",et);
print("alpha=",al);
print("A=",A);
print("T=",T);
print("new_S1=",new_S1(et,om));
print("new_S2=",new_S2(al,et));
print("New S3=",new_S3(al,et,om));
print("S4=",S4(al,T));
print("new_S5=",new_S5(et,om));
print("New S6=",new_S6(al,et,om,A,T));

print("R=",new_R(al,et,om,A,T));
/*
sm=;
print("Sum over zeros came to ",sm);
print("Margin=",1+new_R(al,et,om,A,T)+sm);
f=solve(f=1,4,1+new_R(al,et,om,A,T)+sm+tails(al,et,et/f,om));
f=floor(f*100000);
print("Reducing to eta_0 by a factor of ",f,"/100000");
f=f/100000.0;
print("Tail integrals bounded by ",tails(al,et,et/f,om));
m=-1-new_R(al,et,om,A,T)-sm-tails(al,et,et/f,om);
print("With new et0, Margin=",m);

print("om-et0=",floor((om-et/f)*10000000000),"/10000000000");
print("om+et0=",ceil((om+et/f)*10000000000),"/10000000000");

print(m*exp(floor((om-et/f)*10000000000)/10000000000/2)/ceil((om+et/f)*10000000000)*10000000000);

*/