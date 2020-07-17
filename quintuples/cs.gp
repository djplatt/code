for(r=2,1000,ab=r*r-1;fordiv(ab,a,if(a>=r-1,break);b=ab/a;if(b<a+a,for(c=b+1,1000000,if(issquare(a*c+1)&&issquare(b*c+1),printf("%d %d %d %d\n",r,a,b,c))))));

quit;
