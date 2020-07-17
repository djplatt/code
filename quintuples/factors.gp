for(r=lim1,lim2,f=factor(r*r-1);s=matsize(f);printf("%d %d",r,s[1]);for(r=1,s[1],printf(" %d %d",f[r,1],f[r,2]));printf("\n"));

quit;