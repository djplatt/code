A=30610046000
T=6970346000
et=933831/2^44
al=1153308722614227968

om=727.951332655
/*
al=A^2/om
et=2*A/al
*/

if(om-et<400,print("om-et<400"))
if(4*A/om>al,print("4A/omega>alpha"))
if(al>A^2,print("alpha>A^2"))
if(2*A/al>et,print("2A/alpha>eta"))
if(et>om/2,print("eta>omega/2"))
R1(om,et)=1.5423e-9
R2(al,et,T)=0.08*sqrt(al)*exp(-al*et^2/2)+exp(-T^2/(2*al))*(al/(Pi*T^2)*log(T/(2*Pi))+8*log(T)/T+4*al/T^3)
R3(om,et)=exp((et-om)/2)*log(2*Pi)+3*exp((et-om)/6)
R4(al,om,et,A)=A*log(A)*exp(-A^2/(2*al)+(om+et)/2)*(4*al^(-0.5)+15*et)

print("R1=",R1(om,et))
print("R2=",R2(al,et,T))
print("R3=",R3(om,et))
print("R4=",R4(al,om,et,A))

print("All errors=",R1(om,et)+R2(al,et,T)+R3(om,et)+R4(al,om,et,A))

s=-1.0013360278
print(-1-R1(om,et)-R2(al,et,T)-R3(om,et)-R4(al,om,et,A)-s)

K(y,al)=sqrt(al/(2*Pi))*exp(-al*2*y^2)

tails(om,et,et0,al)=1.3082e-9*K(et0,al)*(et-et0)*(exp((om+et)/2)+exp((om-et0)/2))

print("best ratio for et0=",solve(r=1,10,tails(om,et,et/r,al)+s+1))


