PQa(P,Q,D,root_D,t,u)={\\
   local(a,i,G2,G1,B2,B1,ei,G,B);
   print("In PQa with ",P," ",Q," ",D);\\
   a=floor((P+root_D)/Q);\\
   i=0;\\
   G2=Q;G1=a*Q-P;\\
   B2=0;B1=1;
   while((i==0)||((Q!=1)&&(Q!=-1)),\\
      P=a*Q-P;\\
      Q=(D-P^2)/Q;\\
      ei=(P+root_D)/Q;\\
      if((Q!=1)&&(Q!=-1),\\
         if(ei>1,eib=(P-root_D)/Q;if((eib<0)&&(eib>-1),return([0,0]))));\\
      a=floor(ei);\\
      G=a*G1+G2;\\
      G2=G1;G1=G;\\
      B=a*B1+B2;\\
      B2=B1;B1=B;\\
      i++;\\
      /*print(i," ",P," ",Q," ",a," ",B," ",G);if(i>10,break);*/\\
      );\\
   if(G2^2-D*B2^2>0,return([G2,B2]),if(t!=0,return([G2*t+B2*u*D,B2*t+G2*u]),return([0,0])));\\
}     

PQm1(D,root_D)={\\
   local(P,Q,a,i,G,G1,G2,B,B1,B2);
   P=0;Q=1;\\
   a=floor((P+root_D)/Q);\\
   i=0;\\
   G2=Q;G1=a*Q-P;\\
   B2=0;B1=1;
   while((Q!=1)||(i==0),\\
      P=a*Q-P;\\
      Q=(D-P^2)/Q;\\
      a=floor((P+root_D)/Q);\\
      G=a*G1+G2;\\
      G2=G1;G1=G;\\
      B=a*B1+B2;\\
      B2=B1;B1=B;
      i++;\\
      /*print(i," ",P," ",Q," ",a," ",B," ",G);if(i>10,break);*/\\
      );\\
   if((i%2)==1,return([-1,G2,B2,i]),return([1,G2,B2,i]));\\
}     

cptr=1;
C_LIM=100;
cs=vector(C_LIM);

generate_sols(u,v,x,y,D,ylim,debug)={
   local(x1,y1,t,y2);
   if(debug,print("In generate_sols with ",u," ",v," ",x," ",y," ",D));
   x1=x;y1=y;y2=y;
   while(y2<ylim,
      if(debug,print("x ",x1," y ",y2));
      cs[cptr]=y2;
      cptr++;
      if(cptr==C_LIM,print("Error, cs overflowed.");x=1/0);
      t=x1*u+y1*v*D;
      y1=x1*v+y1*u;
      x1=t;
      t=abs(y1);
      if(y2==t,break,y2=t))}       

solveit(r,a,b,ylim,debug)={
   local(g,a1,b1,D,sfd,N1,N,sN,sD,L1,L2,t);
   if(debug,print("Solving ",a," ",b));
   g=gcd(a,b);\\
   a1=a/g;b1=b/g;\\
   D=a1*b1;sfd=1;\\
   N=(b1-a1)*b1;sN=1;
   sD=sqrt(D);\\
   m1=PQm1(D,sD);\\
   if(m1[1]==-1,
      t=2*m1[2]*m1[3];
      m1[2]=m1[2]*m1[2]+m1[3]*m1[3]*D;
      m1[3]=t);
   L1=0;
   L2=sqrt(N*(m1[2]-1)/(2*D));
   /*printf("L1=%d L2=%d\n",L1,L2);*/
   cptr=1;
   for(y=1,C_LIM,cs[y]=0);
   for(y=0,L2,Sq=N+D*y^2;
      if(issquare(Sq),
         x=floor(sqrt(Sq));
         if(debug,printf("Fund sol found (%d,%d)\n",x,y));
         generate_sols(m1[2],m1[3],x,y,D,ylim,debug);
         generate_sols(m1[2],-m1[3],x,y,D,ylim,debug)));
   scs=vecsort(cs,,8);
   if(debug,print(scs));
   sols=0;
   for(y=0,ylim,t=y*y*D+N;
      if(issquare(t),if(debug,printf("%d ",y));sols++));
   if(debug,print("\n"));
   if(issquare(N),sols--); /* y=0 is a good solution */
   if(sols!=length(scs)-1,
      printf("Mismatch in solutions counts r %d a %d b %d Pell %d BF %d\n",r,a,b,length(scs)-1,sols));
}

for(r=1001,10000,if(r%1000==0,printf("r=%d\n",r));ab=r*r-1;fordiv(ab,a,if(a<r-1,b=ab/a;if(b<=2*a,solveit(r,a,b,1000000,0)))));
