/* elliptic curve over Q of conductor 37 */
r=2
lams=[0.5,0.5]
mus=[0,1]
zeros=[0,5.00195312e+00,6.86914062e+00,8.01367188e+00,9.93164062e+00,1.07753906e+01,1.17558594e+01,1.29589844e+01,1.56035156e+01,1.61933594e+01,1.71425781e+01,1.80644531e+01,1.87871094e+01,1.98144531e+01,2.13222656e+01,2.26191406e+01,2.33300781e+01,2.41699219e+01,2.56582031e+01,2.68144531e+01,2.73378906e+01,2.81894531e+01,2.90292969e+01,2.92832031e+01,3.08964844e+01,3.20410156e+01,3.34394531e+01,3.43652344e+01,3.46347656e+01,3.54628906e+01,3.61621094e+01,3.70839844e+01,3.84667969e+01,3.90019531e+01,3.96035156e+01,4.06503906e+01,4.16542969e+01,4.25761719e+01,4.40332031e+01,4.42363281e+01,4.53925781e+01,4.55644531e+01,4.66542969e+01,4.69589844e+01,4.75761719e+01];
Q=sqrt(37)/Pi
lnG(s)=log(Q)+sum(d=1,length(lams),lams[d]*psi(lams[d]*s+mus[d]))
f(b,h,t)=if(t==0,b/Pi,1-exp(-I*b*t)/(2*Pi*I*t*cosh(h/2*t)))
fhat(b,h,x)=2/Pi*(atan(exp(Pi/h*(x+b)))-atan(exp(Pi/h*(x-b))))
inff(b,h,sig0,k,t)=if(t==0,0,exp(-(lams[k]*sig0/2+mus[k])*t)*(sin(b*lams[k]*t)/(Pi*t))*(1/cosh(h*t*lams[k]/2)-1)/(1-exp(-t)))
wininfint(b,h,sig0,k)=lams[k]*intnum(t=0,1000,inff(b,h,sig0,k,t));
winf1(b,h,sig0)=2*f(b,h,0)*lnG(sig0/2)-sum(k=1,length(lams),wininfint(b,h,sig0,k))
wf(cm,p,m,b,h,sig0)=-cm*sin(b*m*log(p))/(sqrt(p^(m*sig0))*Pi*m*cosh(h*m*log(p)/2))
inff1(b,h,sig0,mu,t)=exp(-sig0/4+mu)*sin(b*t/2)*(1/cosh(h*t/4)-1)/(t*1-exp(-t))
wininfint1(b,h,sig0,mu)=intnum(t=0,1000,inff1(b,h,sig0,mu,t));
winf2(b,h,sig0)=1/Pi*sum(k=1,length(lams),wininfint1(b,h,sig0,mus[k]));
