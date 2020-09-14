/* first two zeros */

z0=14.134725141; z1=21.022039638;

A(t)=0.11+0.29*log(log(t))/log(t)+2.29/log(t)+1/(5*t*log(t));

phi1(t)=1/(t*z0*(z0+t))+1/(2*Pi*t^2)*(dilog(1-(t+z1)/t)+Pi^2/6+log(t)^2/2-log(z1)^2/2-log(2*Pi)*log((t+z1)/z1)+log(z1)*log((t+z1)/t))+A(z1)/t*(2/(z1*(z1+t))*log(z1)+(log(z1/(t+z1))/t^2+1/(t*z1)));

S2=2*phi1(z0)+1/Pi*intnum(t=z1,+oo,log(t/(2*Pi))*phi1(t))+2*A(z1)*(2*phi1(z1)*log(z1)+intnum(t=z1,+oo,phi1(t)/t));

print(S2);


