/*
Lehman's theorem with new s_1

See Saouter & Demichel 2010

*/

/* Assume RH
   T is height to which we have zeros
*/

rh_primepi_low(lx,lplim)=if(lx<=lplim,primepi(exp(lx)),-eint1(-lx)-exp(lx/2)*lx/(8*Pi));
rh_primepi_high(lx,lplim)=if(lx<=lplim,primepi(exp(lx)),-eint1(-lx)+exp(lx/2)*lx/(8*Pi));

no_rh_primepi_low(lx,lplim)=if(lx<=lplim,primepi(exp(lx)),exp(lx)/lx*(1+1/lx+2/(lx*lx)));
no_rh_primepi_high(lx,lplim)=if(lx<=lplim,primepi(exp(lx)),exp(lx)/lx*(1+1/lx+2.334/(lx*lx)));

sd_S1(et,om)=2/(om-et)+10.04/(om-et)^2+log(2)*(om+et)*exp(-(om-et)/2)+2/log(2)*(om+et)*exp(-(om-et)/6);

pi_star_int1(al,et,om,plim,rh)={
   local(lplim);
   lplim=log(plim);
   return(sum(m=2,3,1/m*intnum(u=om-et,om+et,K(u-om,al)*exp(-u/2)*u*no_rh_primepi_high(u/m,lplim))));
}
/*
   local(mlim,res_low,res_high,Kint,lplim);
   lplim=log(plim);
   Kint=intnum(u=om-et,om+et,K(om-u,al));
   if(rh>0,
   res_high=Kint*sum(m=2,floor((om+et)/log(2)),rh_primepi_high((om+et)/m,lplim)/m);
   res_low=Kint*sum(m=2,floor((om-et)/log(2)),rh_primepi_low((om-et)/m,lplim)/m),
   res_high=Kint*sum(m=2,floor((om+et)/log(2)),no_rh_primepi_high((om+et)/m,lplim)/m);
   res_low=Kint*sum(m=2,floor((om-et)/log(2)),no_rh_primepi_low((om-et)/m,lplim)/m));
   
   print("I_1 in [",res_low,",",res_high,"]");
   print("Width of interval =",res_high-res_low);
   return(res_high);
}
*/

my_S1(al,et,om,plim,rh)=pi_star_int1(al,et,om,plim,rh)+log(2)*exp(-(om-et)/2)*(om-et)-1;

S2(al,et)=2*exp(-al*et^2/2)/sqrt(2*Pi*al)/et;

S3(al,et)=0.08*sqrt(al)*exp(-al*et^2/2);

S4(al,T)=exp(-T^2/2/al)*(al/Pi/T^2*log(T/2/Pi)+8*log(T)/T+4*al/T^3);

S5(et,om)=0.05/(om-et);

S6(al,et,om,A,T)=A*log(A)*exp(-A^2/2/al+(om+et)/2)*(4/sqrt(al)+15*et);

R(al,et,om,A,T,plim,rh)=if(rh>0,my_S1(al,et,om,plim,rh)+S2(al,et)+S3(al,et)+S4(al,T)+S5(et,om),my_S1(al,et,om,plim)+S2(al,et)+S3(al,et)+S4(al,T)+S5(et,om)+S6(al,et,om,A,T));

K(y,al)=sqrt(al/2/Pi)*exp(-al*y*y/2);


T=1394846000.0;
om=727.9513330766;
al=2^57;
et=69722785/2251799813685248;
plim=2953652287; /* Dusart works past here */

print("Looking at new region.");
sm=-1.0026225;
/*
print("with RH R=",R(al,et,om,A,T,plim,1));
print("w/o  RH R=",R(al,et,om,A,T,plim,0));
*/

print("Checking Saouter and Demichel");
A=6.85e7;
al=6e12;
T=10946000.0
om=727.961335792;
et=2*A/al;

print("R=",R(al,et,om,A,T,plim,0));
sm=-1.0029061198;
print("Margin=",1+R(al,et,om,A,T,plim,0)+sm);

