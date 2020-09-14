li(lx)=-real(eint1(-lx));

pip(lx)=li(lx)+6.2*lx^1.009*exp(lx-0.83742*sqrt(lx));
pim(lx)=li(lx)-6.2*lx^1.009*exp(lx-0.83742*sqrt(lx));

g(lx)=pip(lx)^2-exp(lx+1)/lx*pim(lx-1);

for(n=55,10000,if(g(n)>=0,printf("%d ",n)));printf("\n");