\p60
nabd=0;
nab=0;
/* given a,b,d and sqrt(b*d+1) find all quadruples {a,b,c,d} with b<c<d */
abd(a,b,d,r)=nabd++;return(abdm(a,b,d,r));

abdp(a,b,d,r)={\\
   local(sols,c,lc,uc,t,s,t1);\\
   return(0);\\
   sols=0;c=0;lc=b+1;uc=d-1;t=1;s=1;\\
   while(1,\\
      t1=t*r+s*d;s=b*t+r*s;t=t1;c=(s*s-1)/b;\\
      if(c>=lc,\\
         if(c<=uc,\\
            /*print("Trying ",a," ",b," ",c," ",d);*/\\
            if(issquare(a*c+1),\\
               sols++;print("Solution:- ",a," ",b," ",c," ",d)),\\
            break)));\\
   return(sols)};

abdm(a,b,d,r)={local(sols,c,lc,uc,t,s,t1);sols=0;c=0;lc=b+1;uc=d-1;t=1;s=-1;while(1,t1=t*r+s*d;s=b*t+r*s;t=t1;c=(s*s-1)/b;if(c>=lc,if(c<=uc,if(issquare(a*c+1),sols++;print("Solution:- ",a," ",b," ",c," ",d)),break)));return(sols)};


/* given a,b and r=sqrt(a*b+1)) find all triples {a,b,d} with b^5<=d<=b^8 */
ab(a,b,r)=nab++;return(abp(a,b,r)+abm(a,b,r));

abp(a,b,r)={local(sols,ld,ud,t,s,d,t1);sols=0;ld=b^5;ud=b^8;t=1;s=1;d=0;while(1,t1=t*r+s*b;s=a*t+r*s;t=t1;d=(s*s-1)/a;if(d>=ld,if(d<=ud,print("a,b,d(+)=",a," ",b," ",d,);sols+=abd(a,b,d,floor(sqrt(b*d+1))),break)));return(sols)};

abm(a,b,r)={local(sols,ld,ud,t,s,d,t1);sols=0;ld=b^5;ud=b^8;t=1;s=-1;d=0;while(1,t1=t*r+s*b;s=a*t+r*s;t=t1;d=(s*s-1)/a;if(d>=ld,if(d<=ud,print("a,b,d(-)=",a," ",b," ",d);sols+=abd(a,b,d,floor(sqrt(b*d+1))),break)));return(sols)};



doit(lim1,lim2)={nabd=0;nab=0;nsols=0;for(r=lim1,lim2,ab1=r*r-1;fordiv(ab1,a,if(a>=r-1,break,b=ab1/a;if(b<a+a,if(b<=1300000000,nsols+=ab(a,ab1/a,r))))));print("There were ",nab," ab pairs.");print("There were ",nabd," abd triples.");print("There were ",nsols," abcd quadruples with r in [",lim1,",",lim2,"]");}

