issq(a,b,r)={return(issquare((r+1)/2));}

compc(a,b,d)={r=floor(sqrt(b*d+1));if(issq(b,d,r),print("Square passed to compc:- ",a," ",b," ",d));c=((d-r)^2-1)/d;\\
/*print("Trying ",a," ",b," ",c," ",d);*/\\
if(issquare(a*c+1),print(a," ",b," ",c," ",d));}

comp_sq(a,b,r)={\\
A=2*a*r;B=2*b*r;R=2*r*r-1;\\
clow=B^5;chi=B^8;\\
if((((r+a)^2-1)%A)!=0,return);\\
x=1;y=1;c=0;\\
while(c<clow,x1=x*r+y*a;y=x*b+y*r;x=x1;c=(x^2-1)/A);\\
while(c<chi,compc(A,B,c);x1=x*r+y*a;y=x*b+y*r;x=x1;c=(x^2-1)/A);\\
x=1;y=-1;c=0;
while(c<clow,x1=x*r+y*a;y=x*b+y*r;x=x1;c=(x^2-1)/A);\\
while(c<chi,compc(A,B,c);x1=x*r+y*a;y=x*b+y*r;x=x1;c=(x^2-1)/A);}

