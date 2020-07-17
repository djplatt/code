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
      print(i," ",P," ",Q," ",a," ",B," ",G);
      );\\
   if((i%2)==1,return([-1,G2,B2,i]),return([1,G2,B2,i]));\\
}     
