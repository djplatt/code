print_facs(q)=
{
  f=factor(q);
  for(n=1,matsize(f)[1],printf("%d ",f[n,1]));
   printf("0\n");
};

do_q(q)=
{
   if((q%2)==1,if((q%4)==1,if(issquarefree(q),return(1),return(0)),return(0)));
   if((q%4)==2,return(0));
   if((q%8)==4,if(issquarefree(q/4),if(((q/4)%4)==3,return(1),return(0)),return(0)));
   if((q%16)==0,return(0));
   q1=q/8;
   if(issquarefree(q1),return(1));
   return(0);
};

print(sum(q=4,1000000000,do_q(q)));

quit;
