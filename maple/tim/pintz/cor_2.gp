e0(lx)=if(lx>=10000,5.2838e-30,if(lx>=9500,4.1547e-29,if(lx>=9000,3.6573e-28,if(lx>=8500,3.1104e-27,if(lx>=8000,2.8158e-26,if(lx>=7500,2.7624e-25,if(lx>=7000,3.1238e-24,if(lx>=6500,3.4824e-23,if(lx>=6000,4.2642e-22,if(lx>=5500,6.0246e-21,if(lx>=5000,8.5638e-20,if(lx>=4500,1.4067e-18,if(lx>=4000,2.3474e-17,if(lx>=3500,4.3021e-16,if(lx>=3000,8.4473e-15,if(lx>=2500,1.8088e-13,if(lx>=2000,4.1103e-12,if(lx>=1500,1.0660e-10,if(lx>=1000,3.7186e-9,if(lx>=900,1.2601e-9,if(lx>=800,1.3065e-9,if(lx>=700,1.3525e-9,if(lx>=600,1.3986e-9,if(lx>=500,1.4451e-9,if(lx>=400,1.4920e-9,if(lx>=300,1.5386e-9,if(lx>=200,1.5859e-9,if(lx>=100,1.6351e-9,if(lx>=90,1.6403e-9,if(lx>=80,1.6495e-9,if(lx>=70,1.6553e-9,if(lx>=60,1.665e-9,if(lx>=50,2.4075e-9,if(lx>=45,1.0412e-8,if(lx>=40,8.6386e-8,1/0)))))))))))))))))))))))))))))))))));

snaplx(lx)=if(lx>10000,10000,if(lx>1000,500*floor(lx/500),if(lx>100,100*floor(lx/100),if(lx>50,50*floor(lx/50),if(lx>40,45,40)))));

f(lx)=6.2*exp(lx)*lx^1.009*exp(-0.83742*sqrt(lx));

/* for lx in [1000,10000] */
E1(lx)={local(slx,res,llx);slx=snaplx(lx);res=0;llx=1500;while(llx<slx,res+=e0(llx-500)*intnum(w=llx-500,llx,exp(w)/w^2);llx+=500);return(res+e0(slx)*intnum(w=slx,lx,exp(w)/w^2)+exp(lx)/lx*e0(slx)+7.6+1/2*intnum(w=log(563),1000,exp(w)/w^3))};

for(lx=1000,10000,if(f(lx+1)<=E1(lx),print(lx)));

/* for lx in [50,100] */
E2(lx)={local(slx,res,llx);slx=snaplx(lx);res=0;llx=60;while(llx<slx,res+=e0(llx-10)*intnum(w=llx-10,llx,exp(w)/w^2);llx+=10);return(res+e0(slx)*intnum(w=slx,lx,exp(w)/w^2)+exp(lx)/lx*e0(slx)+7.6+1/2*intnum(w=log(563),50,exp(w)/w^3))};

for(lx=50,100,if(f(lx+1)<=E2(lx),print(lx)));

E3(lx)={local(slx,res,llx);slx=snaplx(lx);res=0;llx=45;while(llx<slx,res+=e0(llx-5)*intnum(w=llx-5,llx,exp(w)/w^2);llx+=5);return(res+e0(slx)*intnum(w=slx,lx,exp(w)/w^2)+exp(lx)/lx*e0(slx)+7.6+1/2*intnum(w=log(563),40,exp(w)/w^3))};

for(lx=40,50,if(f(lx+1)<=E3(lx),print(lx)));
