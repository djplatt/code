/*
circle * circle
n=100;
a=1;b=0;r1=1;
c=2;d=1;r2=0.5;
forstep(th=0,2*Pi,2*Pi/n,x1=a+r1*cos(th);y1=b+r1*sin(th);print(x1," ",y1));
forstep(th=0,2*Pi,2*Pi/n,x1=c+r2*cos(th);y1=d+r2*sin(th);print(x1," ",y1));
forstep(th=0,2*Pi,2*Pi/n,x1=a+r1*cos(th);y1=b+r1*sin(th);forstep(th1=0,2*Pi,2*Pi/n,x=c+r2*cos(th1);y=d+r2*sin(th1);print(x*x1-y*y1," ",x*y1+y*x1)));
quit
*/

/*

1/rectangle

rectangle in 1st quadrant
*/
n=32;
a=1.0;b=2.0;c=0.5;d=1.0;
/*
draw initial rectangle
*/
forstep(x=a,b,(b-a)/n,print(x," ",c));
forstep(x=a,b,(b-a)/n,print(x," ",d));
forstep(y=c,d,(d-c)/n,print(a," ",y));
forstep(y=c,d,(d-c)/n,print(b," ",y));

/*
draw actual result
*/
forstep(x=a,b,(b-a)/n,z=1.0/(x+c*I);print(real(z)," ",imag(z)));
forstep(x=a,b,(b-a)/n,z=1.0/(x+d*I);print(real(z)," ",imag(z)));
forstep(y=c,d,(d-c)/n,z=1.0/(a+y*I);print(real(z)," ",imag(z)));
forstep(y=c,d,(d-c)/n,z=1.0/(b+y*I);print(real(z)," ",imag(z)));

/*
draw recatangular interval result
*/

w1=a^2+c^2;
w2=b^2+d^2
a=a/w2;
b=b/w1;
temp=-d/w1;
d=-c/w2;
c=temp;
forstep(x=a,b,(b-a)/n,print(x," ",c));
forstep(x=a,b,(b-a)/n,print(x," ",d));
forstep(y=c,d,(d-c)/n,print(a," ",y));
forstep(y=c,d,(d-c)/n,print(b," ",y));

quit
