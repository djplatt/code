sols=0;
ab_count=0;

lim1=1;lim2=100;
for(r=lim1,lim2,ab=r^2-1;fordiv(ab,a,if(a>=r-1,break);b=ab/a;if(b<a+a,ab_count++;solveit4(r,a,b))));

print(sols," solutions found for r in [",lim1,",",lim2,"] comprising ",ab_count," ab pairs");

quit;