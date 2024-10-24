print_facs(q)=
{
  f=factor(q);
  for(n=1,matsize(f)[1],printf("%d ",f[n,1]));
   printf("0\n");
};

do_q(q)=
{
   if((q%2)==1,if((q%4)==1,if(issquarefree(q),printf("%d ",q);print_facs(q)));return;);
   if((q%4)==2,return);
   if((q%8)==4,if(issquarefree(q/4),if(((q/4)%4)==3,printf("%d 4 ",q);print_facs(q/4)));return);
   if((q%16)==0,return);
   q1=q/8;
   if(issquarefree(q1),if((q1%4)==1,printf("%d 8 ",q),printf("%d -8 ",q));print_facs(q1));
   return;
};

for(n=10^9-100,10^9,do_q(n));



quit;

