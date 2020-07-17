allocatemem(1000000000);
/* solve x^2-y^2D=+/-4 , D fundamental
   PQA wont work for D<n^2 (D<16 here)
   so hardwire the ones that break (12)
*/
PQA4(D)={
   local(r,U,V,a0,a,P,Q,P1,Q1,t,P2,Q2);
   if(D==5,return([-4,1,1]));
   if(D==12,return([4,4,1]));
   V=1;a0=floor(sqrt(D));P=a0;Q=1;P1=1;Q1=0;a=a0;
   U=a*V;
   V=(D-U*U)/V;
   while((a<=a0),
      r=P*P-D*Q*Q;
      if(r==4,return([4,P,Q]));
      if(r==-4,return([-4,P,Q]));
      /*if(r==-1,return([-1,2*(P*P+Q*Q*D),4*P*Q]));*/
      if(r==-1,return([-1,2*P,2*Q]));
      if(r==1,return([1,2*P,2*Q]));
      a=floor((a0+U)/V);
      P2=a*P+P1;
      Q2=a*Q+Q1;
      P1=P;
      Q1=Q;
      P=P2;
      Q=Q2;
      U=a*V-U;
      V=(D-U^2)/V;);
return([0,0,0]);}

qr(x,y,d)=log(0.5*(x+y*sqrt(d)));

for(t=3,100,tt=t^2-4;d=quaddisc(tt);issquare(tt/d,&l);f=PQA4(d);print(t," ",d," ",l," ",qfbclassno(d)," ",f[2]," ",f[3]));


/*
print(t," ",d," ",l," ",qfbclassno(d)," ",f[2]," ",f[3]," ");

quit;
*/