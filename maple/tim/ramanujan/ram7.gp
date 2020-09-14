lli(lx)=real(-eint1(-lx));

wh(lx0,lx1,eps)=(1+eps)*(lli(lx1)-lli(lx0))+2*eps*exp(lx0)/lx0;
wl(lx0,lx1,eps)=(1-eps)*(lli(lx1)-lli(lx0))-2*eps*exp(lx1)/lx1;

/* Buthe Gourdon's H interpolating between 100 and 500*/
leps1(lx)=if(lx>=1000,1.94751e-12,if(lx>=500,1.99986e-12,if(lx>=100,1.99986e-12+(500-lx)*(2.45229e-12-1.99986e-12)/400,if(lx>=90,2.52129e-12,if(lx>=75,2.70358e-12,if(lx>=70,2.79233e-12,if(lx>=65,3.57125e-12,if(lx>=60,1.22147e-11,if(lx>=55,1.12494e-10,print("error in leps");break)))))))));

/*
leps(lx)=if(lx>=8000,5.044e-26,if(lx>=7000,5.652e-24,if(lx>=6000,8.528e-22,if(lx>=5000,1.650e-19,if(lx>=4000,4.088e-17,if(lx>=3900,7.020e-17,if(lx>=3800,1.230e-16,if(lx>=3700,2.177e-16,if(lx>=3600,3.810e-16,if(lx>=3500,6.619e-16,if(lx>=3490,7.002e-16,if(lx>=3450,8.783e-16,if(lx>=3400,1.170e-15,if(lx>=100,1.75186e-10,if(lx>=90,1.79331e-10,if(lx>=80,1.6496e-9,if(lx>=70,1.6554e-9,if(lx>=60,1.6666e-9,if(lx>=50,2.4076e-9,print("error in leps");break)))))))))))))))))));
*/

g(lx)={local(ww,eps);eps=leps(lx-1);ww=wh(lx-1,lx,eps);return(ww^2+pim(lx-1)*(pim(lx-1)+2*ww-exp(lx+1)/lx));}

R0=5.573412;
C0=1.0;

rhpip(lx)=lli(lx)+exp(lx/2)/(8*Pi)*lx;
rhpim(lx)=lli(lx)-exp(lx/2)/(8*Pi)*lx;

pim(lx)=if(lx<=55,rhpim(lx),pim(lx-1)+wl(lx-1,lx,leps(lx-1)));
pip(lx)=if(lx<=55,rhpip(lx),pip(lx-1)+wh(lx-1,lx,leps(lx-1)));

gg(lx)=pip(lx)^2-exp(lx+1)/lx*pim(lx-1);
/*

pim(lx)=lli(lx)-7.3e-10*exp(lx)/lx;
pip(lx)=lli(lx)-7.3e+10*exp(lx)/lx;


pim(lx)=lli(lx)-1.001*(sqrt(8/(Pi*R0))*exp(lx)*lx^(-0.75)*exp(-C0*sqrt(lx/R0)));
pip(lx)=lli(lx)+1.001*(sqrt(8/(Pi*R0))*exp(lx)*lx^(-0.75)*exp(-C0*sqrt(lx/R0)));

*/


read("../pintz/pintz_nocalc.gp");
if(type(lexps)!="t_VEC",lexps=vector(4000);for(n=30,79,lexps[n*50]=chk(n*50);for(m=n*50+1,n*50+49,lexps[m]=lexps[n*50]));lexps[4000]=chk(4000);for(n=3395,3403,lexps[n]=chk(n)));


leps(lx)=if(lx>=1500,lexps[lx],if(lx>=1000,1.43770e-10,if(lx>=500,1.47067e-10,if(lx>=100,1.75185e-10,if(lx>=90,1.79330e-10,if(lx>=80,1.84848e-10,if(lx>=70,1.91910e-10,if(lx>=65,1.96865e-10,if(lx>=60,2.08162e-10,if(lx>=55,2.88434e-10,if(lx>=log(59),lx^2/8/Pi/exp(lx/2),print("error in leps");break)))))))))));

pims=vector(4000);
pims[55]=pim(55);
pims[54]=pim(54);
for(n=56,4000,pims[n]=pims[n-1]+wl(n-1,n,leps(n-2)));
pips=vector(4000);
pips[55]=pip(55);
pips[54]=pip(54);
for(n=56,4000,pips[n]=pips[n-1]+wh(n-1,n,leps(n-2)));
ggg(lx)=pips[lx]^2-pims[lx-1]*exp(lx+1)/lx;

for(n=56,4000,if(ggg(n)>=0,printf("%d ",n)));printf("\n");
