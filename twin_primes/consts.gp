x0=4e18;
pi2x0=3023463123235320;
B0=1.840518;


al=2/5;
Hal=950.05;
cal=1.0503; /* old c if one considers x in (0,1]. */
cal=0.6877; /* this happens as we approch x=6 */

/*
al=1/3;
Hal=251.02;
cal=1.641;
cal=;
*/
/* First term from H.cpp, second from H_tail_estimate.mw */

/*
al=7/16;
Hal=3412.331;
*/

rho=sqrt(1+2/3*sqrt(6/5));
A6=9.27436-2*log(rho);
A7=-5.6646+log(rho)^2-9.2744*log(rho);
A9=24.09391*sqrt(rho);
C=1.320324;
A8=16*C*cal*Hal*rho^(al/2);

myroundup(x,n)=ceil(x*n)/n;
myrounddown(x,n)=floor(x*n)/n;
A6=myrounddown(A6,100000);
A7=myrounddown(A7,100000);
A8=myroundup(A8,100000);
A9=myroundup(A9,10000000);

printf("A6=%e\n",A6);
printf("A7=%e\n",A7);
printf("A8=%e\n",A8);
printf("A9=%e\n",A9);




F(x,lx,sx)=A6+A7/lx-A8/x^(al/2)/lx-A9/(sx*lx);
print("F(exp(100))=",F(exp(100),100,exp(50)));

LL=690;

/* Table 3.4 of Klyve */
D(L)=if(L>=690,8.45,if(L>=396,8.44,if(L>=278,8.43,if(L>=214,8.42,if(L>=174,8.41,if(L>=147,8.4,if(L>=127,8.39,if(L>=100,8.37,if(L>=82,8.35,if(L>=60,8.3,if(L>=48,8.2,if(L>=44,8.1,print("log x too small in D()";break)))))))))))));



/* Note Klyve has 16*C */
pi2_byx(lx)=8*C/(lx*(D(lx)+lx));
/*
print(B0-2*pi2x0/x0+2*intnum(x=x0,exp(690),pi2_byx(log(x))/x)+16*C/690.);
*/

/* 3.20 of Riesel & Vaughan */
pi2(x)={local(lx,sx);lx=log(x);sx=sqrt(x);return(2*sx+8*C*x/(lx*(lx+F(x,lx,sx))));}



N=1000;
II=2*intnum(t=x0,exp(50),pi2(t)/t^2)+2*sum(n=5,N-1,intnum(t=exp(n*10),exp((n+1)*10),pi2(t)/t^2));
print("x0=",x0," Integral = ",II);
print(B0-2*pi2x0/x0+II+16*C/(N*10));
quit;

/*
x0=10^19;
B0=1.84181;
pi2x0=7.2376e15;
II=2*intnum(t=x0,exp(50),pi2(t)/t^2)+2*sum(n=5,N,intnum(t=exp(n*10),exp((n+1)*10),pi2(t)/t^2));
printf("x0=%e Integral = %e\n",x0,II);
print(B0-2*pi2x0/x0+II+16*C/(N*10));

x0=10^20;
B0=1.84482;
pi2x0=6.5155e16;
II=2*intnum(t=x0,exp(50),pi2(t)/t^2)+2*sum(n=5,N,intnum(t=exp(n*10),exp((n+1)*10),pi2(t)/t^2));
printf("x0=%e Integral = %e\n",x0,II);
print(B0-2*pi2x0/x0+II+16*C/(N*10));

x0=10^80;
B0=1.8878;
pi2x0=3.9341e75;
II=2*intnum(t=x0,exp(190),pi2(t)/t^2)+2*sum(n=19,N,intnum(t=exp(n*10),exp((n+1)*10),pi2(t)/t^2));
printf("x0=%e Integral = %e\n",x0,II);
print(B0-2*pi2x0/x0+II+16*C/(N*10));


N=1000;
x0=4e18;
pi2x0=3023463123235320;
B0=1.840518;
A6=9.27436;
A7=0;
A9=0;
A8=0;
II=2*intnum(t=x0,exp(50),pi2(t)/t^2)+2*sum(n=5,N,intnum(t=exp(n*10),exp((n+1)*10),pi2(t)/t^2));
print("Integral = ",II);
print(B0-2*pi2x0/x0+II+16*C/(N*10));
*/
