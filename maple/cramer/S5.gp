A(t)=0.11+0.29*log(log(t))/log(t)+2.29/log(t)+1/(5*t*log(t));

phi5(t)=1/(8*Pi*t)*log((t+2)/t)*log(t*(t+2)/(4*Pi*Pi))+A(t)/t*(log(t)/t+1/(t*(t+2)));

S2=2*(1/(2*Pi)*intnum(t=26000,+oo,log(t/(2*Pi))*phi5(t))+A(26000)*(2*phi5(26000)*log(26000)+intnum(t=26000,+oo,phi5(t)/t)));

print(S2);
