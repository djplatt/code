
/*
The original e0 function to use when log(x) is less than 35
*/
e00(x)={local(R,X);R=6.455;X=sqrt(log(x)/R);sqrt(8*X/(17*Pi))*exp(-X)}

/*
e00(x)=print("This should never get called. ",log(x));

e0(x)={local(R,X,lx);lx=log(x);if(lx>=100,1.75185e-10,if(lx>=90,1.79330e-10,if(lx>=80,1.84848e-10,if(lx>=70,1.9191e-10,if(lx>=65,1.96865e-10,if(lx>=55,2.88434e-10,if(lx>=50,1.16465e-9,if(lx>=46,3.46e-9,if(lx>=45,4.82e-9,if(lx>=44,6.62e-9,if(lx>=43,9.58e-9,if(lx>=42,1.45e-8,if(lx>=41,2.26e-8,if(lx>=40,3.63e-8,if(lx>=39,5.26e-8,if(lx>38,7.94e-8,if(lx>=37,1.24e-7,if(lx>=36,1.97e-7,if(lx>=34.99999999,3.08e-7,e00(x))))))))))))))))))));}
*/
e0(x)={local(R,X,lx);lx=log(x);if(lx>=100,1.75185e-10,if(lx>=90,1.79330e-10,if(lx>=80,1.84848e-10,if(lx>=70,1.9191e-10,if(lx>=65,1.96865e-10,if(lx>=55,2.88434e-10,if(lx>=50,1.16465e-9,if(lx>=45,4.82e-9,if(lx>=40,3.63e-8,if(lx>=34.99999999,3.08e-7,e00(x)))))))))));}

e0(x)=if(log(x)>=35,3.08e-7,e00(x))

lx0=35;
x0=exp(lx0);
I_hi=46690548661120-x0/lx0*(1-e0(x0));
I_low=46690548661120-x0/lx0*(1+e0(x0));

pip(x)=if(x<=1.39e17,x/log(x)+intnum(t=2,x,1/log(t)^2),x*(1+e0(x))/log(x)+I_hi+intnum(t=x0,x,(1+e0(t))/log(t)^2));
pip(x)=x*(1+e0(x))/log(x)+I_hi+intnum(t=log(x0),log(x),exp(t)*(1+e0(exp(t)))/t^2);
pim(x)=x*(1-e0(x))/log(x)+I_low+intnum(t=log(x0),log(x),exp(t)*(1-e0(exp(t)))/t^2);

g(x)=pip(x)^2-exp(1)*x/log(x)*pim(x/exp(1));


print(36," ",g(exp(36)));
print(9658," ",g(exp(9658)));


