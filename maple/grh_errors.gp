allocatemem(2^30);
A=64/5;
QT=1e8;
turing_region=30;
T_FAC=7.5;
/* odd characters */
E_FAC=7.03; /*7.03; /* 4.0 for even */
E_FACe=4.0;

/*
X(x)=Pi*delta*exp(-delta)*exp(x)^2/q;

f_hat_twiddle_err_odd(n)={local(w1,w2,x1,x2);w1=2*Pi*(A+n/B);w2=2*Pi*(A-n/B);x1=X(w1);if(x1<=1,print("X(w1) too small");return(0));x2=X(w2);if(x2<=1,print("X(w2) too small");return(0));return((4*exp(3*w1/2-X(w1))*(1+1/(2*X(w1)))^1.5+exp(3*w2/2-X(w2))*(1+1/(2*X(w2)))^1.5)/(q^0.75*delta^0.5*(1-exp(-Pi*A))));}
*/

Eo(t)=zeta(9/8)*Pi^-0.75*exp(real(lngamma(3/4+I*t/2))+Pi*t*eta_/4)*(q/(2*Pi)*abs(1.5+t))^(5/16);

Betao(t)=Pi/4-1.5*atan(1/2/abs(t))-4/(Pi*Pi*abs(t^2-9/4));

f_twiddle_err_odd(m)=Eo(m/A+B)/(1-exp(-B*(Betao(m/A+B)-Pi*eta_/4)))+Eo(m/A-B)/(1-exp(-B*(Betao(m/A-B)+Pi*eta_/4)));

Ee(t)=zeta(9/8)*Pi^-0.25*exp(real(lngamma(1/4+I*t/2))+Pi*t*etae/4)*(q/(2*Pi)*abs(1.5+t))^(5/16);

Betae(t)=Pi/4-0.5*atan(1/2/abs(t))-4/(Pi*Pi*abs(t^2-1/4));

f_twiddle_err_even(m)=Ee(m/A+B)/(1-exp(-B*(Betae(m/A+B)-Pi*eta_/4)))+Ee(m/A-B)/(1-exp(-B*(Betae(m/A-B)+Pi*eta_/4)));

h=7/32;
Ni=20;

G(t0,n)=(3/2+t0+(Ni+n)/A)^(9/16)*exp(-(Ni+n)^2/(2*A^2*h^2))/(Pi*(Ni+n));

int_err()={sqrt(Pi)*zeta(9/8)*exp(1/6)*2^(5/4)*(q/(2*Pi))^(5/16)*G(T0,0)/(1-G(T0,1)/G(T0,0))};

set_N(T0)={local(p2);p2=1;while(p2<(T0*T_FAC)*A,p2=p2+p2);return(p2)};

qs=[3,11,100,1000,10000,87381,87382];

for(qi=3,10000,q=qi;T0=ceil(QT/(10.0*q))*10.0+turing_region;N=set_N(T0);B=N/A;eta_=1-E_FAC/(QT/q+turing_region);etae=1-E_FACe/(QT/q+turing_region);print("q=",q," N=2^",round(log(N)/log(2))," fte(even)=",f_twiddle_err_even(0)));

quit;

/*
print("q=",q," inter error =",int_err()," fte(even)=",f_twiddle_err_even(0)," fte(odd)=",f_twiddle_err_odd(0)));
*/