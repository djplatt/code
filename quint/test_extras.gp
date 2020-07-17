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

solveit(r,a,b)={\\
   g=gcd(a,b);\\
   a1=a/g;b1=b/g;\\
   D=/*core(a1*b1)*/a1*b1;sfd=floor(sqrt(a1*b1/D));\\
   N=(b1-a1)*b1;sD=sqrt(D);\\
   m1=PQm1(D,sD);\\
   /*print("m1=",m1);*/\\
   if(m1[1]==-1,t=m1[2];u=m1[3],t=0;u=0);\\
   fordiv(N,f,if(issquare(f),F=floor(sqrt(f));m=N/f;if(m%2==1,m2=(m-1)/2;m1=-m2,m2=m/2;m1=-m2+1);for(z=m1,m2,if(((z^2%m)==(D%m)),V=PQa(z,m,D,sD,t,u);V[1]*=F;V[1]/=b1;V[2]*=F;V[2]/=sfd;V[1]=abs(V[1]);V[2]=abs(V[2]);if((V[1]!=0)&&(V[1]!=1)&&(type(V[1])=="t_INT")&&(type(V[2])=="t_INT"),print(r," ",a," ",b," ",t," ",u," ",V," ",b*V[1]^2-a*V[2]^2," ",b-a))))));\\
}

solveit1(r,a,b)={\\
   a1=a;b1=b;\\
   D=a1*b1;sfd=1;\\
   N=(b1-a1)*b1;sD=sqrt(D);\\
   m1=PQm1(D,sD);\\
   /*print("m1=",m1);*/\\
   if(m1[1]==-1,t=m1[2];u=m1[3],t=0;u=0);\\
   fordiv(N,f,if(issquare(f),F=floor(sqrt(f));m=N/f;if(m%2==1,m2=(m-1)/2;m1=-m2,m2=m/2;m1=-m2+1);for(z=m1,m2,if(((z^2%m)==(D%m)),V=PQa(z,m,D,sD,t,u);V[1]*=F;V[1]/=b1;V[2]*=F;V[2]/=sfd;V[1]=abs(V[1]);V[2]=abs(V[2]);if((V[1]!=0)&&(V[1]!=1)&&(V[1]!=r-a)&&(type(V[1])=="t_INT")&&(type(V[2])=="t_INT"),print(r," ",a," ",b," ",t," ",u," ",V," ",b*V[1]^2-a*V[2]^2," ",b-a))))));\\
}

solveit2(r,a,b)={\\
   g=gcd(a,b);\\
   a1=a/g;b1=b/g;\\
   D=core(a1*b1);sfd=floor(sqrt(a1*b1/D));\\
   N1=(b1-a1)*b1;\\
   N=core(N1);\\
   sN=floor(sqrt(N1/N));\\
   sD=sqrt(D);\\
   /*print("Solving ",sN,"(X^2-",D,"Y^2=",N,")");
   print("X=",b1/sN,"x Y=",sfd/sN,"y");*/
   m1=PQm1(D,sD);\\
   /*print("m1=",m1);*/\\
   if(m1[1]==-1,t=m1[2];u=m1[3],t=0;u=0);\\
   fordiv(N,f,if(issquare(f),F=floor(sqrt(f));m=N/f;if(m%2==1,m2=(m-1)/2;m1=-m2,m2=m/2;m1=-m2+1);for(z=m1,m2,if(((z^2%m)==(D%m)),V=PQa(z,m,D,sD,t,u);print("PQa returned ",V);V[1]*=F*sN;V[1]/=b1;V[2]*=F*sN;V[2]/=sfd;V[1]=abs(V[1]);V[2]=abs(V[2]);if((V[1]!=0)&&(V[1]!=1)&&(V[1]!=r-a)&&(type(V[1])=="t_INT")&&(type(V[2])=="t_INT"),print(r," ",a," ",b," ",t," ",u," ",V," ",b*V[1]^2-a*V[2]^2," ",b-a))))));\\
}

solveit3(r,a,b)={\\
   local(D,N,t,L1,count);\\
   D=a*b;N=b*(b-a);\\
   t=r;\\
   L1=sqrt(N*(t-1)/D);\\
   /*print("Solving ",r," ",a," ",b," ",L1);*/\\
   count=0;\\
   for(y=0,L1,Sq=N+D*y^2;if(issquare(Sq),if(floor(sqrt(Sq))%b==0,count++)));\\
   if(count>1,\\
      print(r," ",a," ",b," had ",count," fundamental solutions.");return(1)\\
      ,return(0));\\
}   

bf(N,D,sD)={
   m1=PQm1(D,sD);\\
   if(m1[1]==-1,print(b,"x^2-",a,"y^2=-1 has solution.");return);
   L1=sqrt(N*(m1[2]-1)/D);
   for(y=0,L1,Sq=N+D*y^2;if(issquare(Sq),x=floor(sqrt(Sq));print("try (",m1[2],"+",m1[3],"(",D,")^1/2)^n * (",x,"+",y,"(",D,")^1/2)")));
}

cptr=1;dptr=1;
C_LIM=100;D_LIM=100;
cs=vector(C_LIM);
ds=vector(D_LIM);
generate_cs_ds(u,v,x0,y0,D,a,b,b1)=
{
   local(b2,c,lwb,upb,x,y);
   /*print("In generate_cs_ds with ",u," ",v," ",x0," ",y0," ",a," ",b," ",b1," ",D);*/
   x=x0;y=y0;
   lwb=b^5;upb=b^8;
   b2=b1*b1;
   c=(x*x/b2-1)/a;
   while(c<lwb,
      if((type(c)=="t_INT")&&(c>b),
         cs[cptr]=[c];
         cptr++);
      t=x*u+y*v*D;
      y=x*v+y*u;
      x=t;
      c=(x*x/b2-1)/a);
   /*for(i=1,cptr-1,print("c found=",cs[i]));*/
   d=c;
   while(d<=upb,
      if(type(d)=="t_INT",
         ds[dptr]=[d];
         dptr++);
      t=x*u+y*v*D;
      y=x*v+y*u;
      x=t;
      d=(x*x/b2-1)/a);
}

generate_cs_ds_extra(u,v,x0,y0,D,a,b,b1)=
{
   local(b2,c,lwb,upb,x,y);
   /*print("In generate_cs_ds with ",u," ",v," ",x0," ",y0," ",a," ",b," ",b1," ",D);*/
   x=x0;y=y0;
   lwb=b^5;upb=b^8;
   b2=b1*b1;
   c=(x*x/b2-1)/a;
   while(c<lwb,
      if(c>b,
         cs[cptr]=[c];
         cptr++);
      t=x*u+y*v*D;
      y=x*v+y*u;
      x=t;
      c=(x*x/b2-1)/a);
   d=c;
   while(d<=upb,
      ds[dptr]=[d];
      dptr++;
      t=x*u+y*v*D;
      y=x*v+y*u;
      x=t;
      d=(x*x/b2-1)/a);
}

sols=0;
/* test every c>b against every d */

test_ds(r,a,b)={
   local(c,d,lastc,lastd,scs,sds);
   scs=vecsort(cs,1,8);
   sds=vecsort(ds,1,8);
   for(i=2,length(scs),
      c=scs[i][1];
      if(type(c)=="t_INT",
         for(j=2,length(sds),
            d=sds[j][1];
            if(type(d)=="t_INT",
               if(issquare(c*d+1),
                  printf("Solution %d %d %d %d %d.\n",r,a,b,c,d);
                  sols++;
                  )))));
}
         

solveit4(r,a,b)={
   local(g,a1,b1,D,sfd,N1,N,sN,sD,L1,L2);
   for(i=1,C_LIM,cs[i]=[0]);
   for(i=1,D_LIM,ds[i]=[0]);
   cptr=1;dptr=1;
   g=gcd(a,b);\\
   a1=a/g;b1=b/g;\\
   D=a1*b1;sfd=1;\\
   N=(b1-a1)*b1;sN=1;
   sD=sqrt(D);\\
   m1=PQm1(D,sD);\\
   if(m1[1]==-1,
      /*print("m1=",m1);*/
      L1=ceil(sqrt(N/D));
      L2=floor(sqrt(N*(m1[2]+1)/(2*D)));
      /*printf("L1=%d L2=%d\n",L1,L2);*/
      for(y=L1,L2,Sq=D*y^2-N;
         if(issquare(Sq),
            x=floor(sqrt(Sq));
            /*printf("Fund sol found (%d,%d)\n",x,y);*/
            generate_cs_ds(m1[2],m1[3],x,y,D,a,b,b1);
            generate_cs_ds(m1[2],-m1[3],x,y,D,a,b,b1)));
      t=m1[2]^2+m1[3]^2*D;
      m1[3]*=2*m1[2];
      m1[2]=t;
      /*for(x=1,m1[2],for(y=1,m1[3],if((x^2-D*y^2==1)&&((x<m1[2])||(y<m1[3])),
         printf("Missed solution %d %d %d %d\n",a,b,x,y))));*/);
   L1=0;
   L2=sqrt(N*(m1[2]-1)/(2*D));
   /*printf("L1=%d L2=%d\n",L1,L2);*/
   for(y=0,L2,Sq=N+D*y^2;
      if(issquare(Sq),
         x=floor(sqrt(Sq));
         /*printf("Fund sol found (%d,%d)\n",x,y);*/
         generate_cs_ds(m1[2],m1[3],x,y,D,a,b,b1);
         generate_cs_ds(m1[2],-m1[3],x,y,D,a,b,b1)));
   /*for(i=1,cptr-1,if(cs[i]<=1000000,printf("%d %d %d %d\n",r,a,b,cs[i])));*/
   test_ds(r,a,b);
}


test_extra(a,b,u,v)={
   local(g,a1,b1,D,sfd,N1,N,sN,sD,L1,L2);
   for(i=1,C_LIM,cs[i]=[0]);
   for(i=1,D_LIM,ds[i]=[0]);
   cptr=1;dptr=1;
   g=gcd(a,b);\\
   a1=a/g;b1=b/g;\\
   D=a1*b1;sfd=1;\\
   N=(b1-a1)*b1;sN=1;
   sD=sqrt(D);\\
   m1=PQm1(D,sD);\\
   if(m1[1]==-1,
      print("m1=",m1);
      t=m1[2]^2+m1[3]^2*D;
      m1[3]*=2*m1[2];
      m1[2]=t);
   generate_cs_ds_extra(m1[2],m1[3],u,v,D,a,b,b1);
   generate_cs_ds_extra(m1[2],-m1[3],u,v,D,a,b,b1);   
   test_ds(r,a,b);
}

test_extras(Y)=for(ptr=1,length(Y)/4,a=Y[ptr*4-3];b=Y[ptr*4-2];x=Y[ptr*4-1];y=Y[ptr*4];test_extra(a,b,x,y));

Y=readvec("~/data/quint/extras.log");

printf("Processing %d extra solutions.\n",length(Y)/4);

test_extras(Y);

printf("We had %d solutions.\n",sols);

quit;