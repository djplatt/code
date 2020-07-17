setX(t1,r,mus)={local(X1);X=floor(sqrt(t1^2-(5/2+mus[1])^2));for(j=2,r,X1=floor(sqrt(t1^2-(5/2+mus[j])^2));if(X1<X,X=X1))}
Q(s,N,r,mus)=N*prod(j=1,r,(s+mus[j])/(2*Pi));

Sint(t1,t2,N,r,mus,X)=(1/4*log(abs(Q(3/2+I*t2,N,r,mus)))+(log(2)-1/2)*log(abs(Q(3/2+I*t1,N,r,mus)))+5.65056*r+r/sqrt(2)/(X-5))/Pi;

phi(t,N,r,mus)=(t*log(N)/2-log(Pi)*r*t/2+sum(j=1,r,imag(lngamma((1/2+I*t+mus[j])/2))))/Pi;

phiint(t1,t2,N,r,mus)=intnum(t=t1,t2,phi(t,N,r,mus));

