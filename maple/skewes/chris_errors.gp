K(y,alpha)=sqrt(alpha/(2*Pi))*exp(-alpha*y^2/2);
li(x)=-real(eint1(-log(x)));
R1(om,et,al,T)=1.812*(om+et)/exp((om-et)/6);
R2(om,et,al,T)=0.024*exp(-(om-et)/4)*(1+4/(om-et))+exp(1/32/al-om/4)*(1.301+0.04*al)+2*exp(1/32/al-om/4)/(al*et-1/4)*(log((4*al*et-1)/2/Pi)^2+log(4*al*et-1)+0.9321);
R3(om,et,al,T)=al*exp(-T^2/2/al)/2/Pi*log(T/2/Pi)*(1/T^3+2/T^2);
R4(om,et,al,T)=K(et,al)/al/et*(log(al*et/Pi)^2/Pi+4*log(2*al*et)+4.52)
R5(om,et,al,T)=2*al*exp(-T^2/2/al)*log(T/2/Pi)/(2*Pi*om*T^3)+K(et,al)/(al*om*et)*(0.047+1/al/et)+0.019/sqrt(al)/om^2;
R6(om,et,al,T)=2.92e-3/(om-et)^2;

et=1.3e-9;
om=727.9513329826;
T=3.0610046e10;
al=4.5e19;
for(n=0,5,om=(7279513329824+n)/10^10;\
II=0.5*intnum(u=om-et,om+et,K(u-om,al)*u*exp(-u/2)*li(exp(u/2)));\
r1=R1(om,et,al,T);\
r2=R2(om,et,al,T);\
r3=R3(om,et,al,T);\
r4=R4(om,et,al,T);\
r5=R5(om,et,al,T);\
r6=R6(om,et,al,T);\
print("alpha=",al);\
print("R1=",r1);\
print("R2=",r2);\
print("R3=",r3);\
print("R4=",r4);\
print("R5=",r5);\
print("R6=",r6);\
print("int=",II);\
print("Tot="r1+r2+r3+r4+r5+r6+II);print(""););
