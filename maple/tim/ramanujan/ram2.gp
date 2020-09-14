
/*
The original e0 function to use when log(x) is less than 35
*/
e00(x)={local(R,X);R=6.455;X=sqrt(log(x)/R);sqrt(8*X/(17*Pi))*exp(-X)}

/*
A better one for large x
*/
e01(x)={local(X);X=sqrt(log(x)/5.573412);sqrt(8*X/Pi)*exp(-X)}

/*
A revised one using Buethe/FK bounds
*/
e0(x)={local(R,X,lx);lx=log(x);if(lx>=100,1.75185e-10,if(lx>=90,1.79330e-10,if(lx>=80,1.84848e-10,if(lx>=70,1.9191e-10,if(lx>=65,1.96865e-10,if(lx>=55,2.88434e-10,if(lx>=50,1.16465e-9,if(lx>=45,4.82e-9,if(lx>=44,6.62e-9,if(lx>=43,9.58e-9,if(lx>=42,1.45e-8,if(lx>=41,2.26e-8,if(lx>=40,3.63e-8,if(lx>=39,5.26e-8,e00(x)))))))))))))))};





e0(x)={local(R,X,lx);lx=log(x);if(lx>=3355,e01(x),if(lx>=100,1.75185e-10,if(lx>=90,1.79330e-10,if(lx>=80,1.84848e-10,if(lx>=70,1.9191e-10,if(lx>=65,1.96865e-10,if(lx>=55,2.88434e-10,if(lx>=50,1.16465e-9,if(lx>=46,3.46e-9,if(lx>=45,4.82e-9,if(lx>=44,6.62e-9,if(lx>=43,9.58e-9,if(lx>=42,1.45e-8,if(lx>=41,2.26e-8,if(lx>=40,3.63e-8,if(lx>=39,5.26e-8,if(lx>38,7.94e-8,if(lx>=37,1.24e-7,if(lx>=36,1.97e-7,if(lx>=35,3.08e-7,if(lx>=30,2.811e-6,e00(x))))))))))))))))))))))}



pip(x)=x*(1+e0(x))/log(x)+intnum(t=log(2),log(x),exp(t)*(1+e0(exp(t)))/t^2);
pim(x)=x*(1-e0(x))/log(x)+intnum(t=log(2),log(x),exp(t)*(1-e0(exp(t)))/t^2);
g(x)=pip(x)^2-exp(1)*x/log(x)*pim(x/exp(1));


