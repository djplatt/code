\p400
allocatemem(2^30);
MAX_MU=8;
MAX_MU2=101;
MAX_R=9;
h=4;
trunc(x)=round(x*2^20); /* +/- 1 added in buthe.h */
printf("{");
for(r=2,2,for(kk=0,2*MAX_MU2,k=kk/2.0;b=64.0/r;S=intnum(x=0,10000,2*exp(-(0.5+k)*x)/(1-exp(-2*x))*(b/Pi-sin(b*x)/(Pi*x*cosh(1/2*h*x))));printf("%d,",trunc(S)));printf("\n");)
for(r=3,MAX_R,for(kk=0,2*MAX_MU,k=kk/2.0;b=64.0/r;S=intnum(x=0,10000,2*exp(-(0.5+k)*x)/(1-exp(-2*x))*(b/Pi-sin(b*x)/(Pi*x*cosh(1/2*h*x))));printf("%d,",trunc(S)));for(kk=2*MAX_MU+1,2*MAX_MU2,printf("0,"));printf("\n");)
printf("}\n");
quit;
