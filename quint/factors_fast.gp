for(r=lim1,lim2,f1=factor(r-1);s1=matsize(f1);f2=factor(r+1);s2=matsize(f2);printf("%d %d",r,s1[1]+s2[1]);r1=1;r2=1;while((r1<=s1[1])&&(r2<=s2[1]),if(f1[r1,1]<=f2[r2,1],printf(" %d %d",f1[r1,1],f1[r1,2]);r1++,printf(" %d %d",f2[r2,1],f2[r2,2]);r2++));if(r1<=s1[1],for(p=r1,s1[1],printf(" %d %d",f1[p,1],f1[p,2])),for(p=r2,s2[1],printf(" %d %d",f2[p,1],f2[p,2])));printf("\n"));

quit;
