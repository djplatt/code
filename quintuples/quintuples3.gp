lim1=100000;

/* given a,b,c,d and sqrt(b*d+1) find all triples {a,b,c,d} with b<c<d */
abd(a,b,d,r)=return(abdp(a,b,d,r)+abdm(a,b,d,r));

abdp(a,b,d,r)={sols=0;c=0;lc=b+1;uc=d-1;t=1;s=1;while(1,t1=t*r+s*d;s=b*t+r*s;t=t1;c=(s*s-1)/b;if(c>lc,if(c<=uc,if(issquare(a*c+1),sols++;print(a," ",b," ",c," ",d)),break)));return(sols)};

abdm(a,b,d,r)={sols=0;c=0;lc=b+1;uc=d-1;t=1;s=-1;while(1,t1=t*r+s*d;s=b*t+r*s;t=t1;c=(s*s-1)/b;if(c>lc,if(c<=uc,if(issquare(a*c+1),sols++;print(a," ",b," ",c," ",d)),break)));return(sols)};


/* given a,b and sqrt(a*b+1)) find all triples {a,b,d} with b^6<=d<=b^7.7 */
ab(a,b,r)=return(abp(a,b,r)+abm(a,b,r));

abp(a,b,r)={sols=0;ld=b^6;ud=floor(b^7.7);t=1;s=1;d=0;while(1,t1=t*r+s*b;s=a*t+r*s;t=t1;d=(s*s-1)/a;if(d>=ld,if(d<=ud,print(a," ",b," ",d);sols+=abd(a,b,d,floor(sqrt(b*d+1))),break)));return(sols)};

abm(a,b,r)={sols=0;ld=b^6;ud=floor(b^7.7);t=1;s=-1;d=0;while(1,t1=t*r+s*b;s=a*t+r*s;t=t1;d=(s*s-1)/a;print("s=",s," ","d=",d);if(d>=ld,if(d<=ud,print(a," ",b," ",d);sols+=abd(a,b,d,floor(sqrt(b*d+1))),break)));return(sols)};



doit(lim1,lim2,lim3)={nsols=0;for(a=lim1,lim2,for(b=a+3,lim3,ab1=a*b+1;if(issquare(ab1),nsols+=ab(a,b,floor(sqrt(ab1))))));print("We found ",nsols," solutions with a<b<=",lim1);}

doit(1,1,63);

quit;
