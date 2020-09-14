th1(t1)=(0.1575*log(4*t1)+3.2225)/log(t1);

E1(t1,al)=2*sqrt(2*Pi)*(sqrt(t1^2+1/4)/al)^(1/2/al-1/2)*exp(al/6/sqrt(t1^2+1/4))*exp(-Pi*t1/2/al)*(al/Pi^2*log(2*t1*exp(1)/Pi)+th1(t1)*log(t1^2*exp(1)));
