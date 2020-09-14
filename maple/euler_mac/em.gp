f1(t)=t^(-al);

EM1(x)=1/(1-al)*(x^(1-al)-1)+1/2*x^(-al)-al/12*x^(-al-1)+1/2+al/12-al*(al+1)*(al+2)/6*intnum(t=1,x,(frac(t)^3-1.5*frac(t)^2+0.5*frac(t))*t^(-al-3));

EM2(x)=x^(1-al)*log(x)/(1-al)+1/(1-al)^2*(1-x^(1-al))+0.5*x^(-al)*log(x)-1/12-1/6*intnum(t=1,x,(al*(al+1)*(al+2)*log(t)-(3*al^2+6*al+2))*(frac(t)^3-1.5*frac(t)^2+0.5*frac(t))*t^(-al-3));