/* Pari code to support version 23 */
R=5.573412;

C1s=[5.277,6.918,8.975,9.441,9.926,10.431,10.955,11.499,12.063,12.646,13.250,13.872,14.513,15.173,15.850,16.544,17.253];

C2s=[4.403,3.997,3.588,3.514,3.430,3.346,3.262,3.186,3.102,3.017,2.941,2.856,2.772,2.694,2.609,2.532,2.446];

C(sig,Cs)={sig1=sig*100;if(sig1==75,return(Cs[1]));if(sig1==80,return(Cs[2]));if(floor(sig1)!=sig1,return(0));if(sig1<85,return(0));if(sig1>99,return(0));return(Cs[floor(sig1)-85+3])};

C1(sig)=C(sig,C1s);
C2(sig)=C(sig,C2s);

k0(sig,lx)=1.0/(2.002*C1(sig)*exp(((10-16*sig)/3)*sqrt(lx/R))*(2*sqrt(lx/R))^(5-2*sig));

C3(sig,lx,k)=2*exp(-2*sqrt(lx/R))*lx^2*k;

C4(sig,lx,k)=exp(lx*(sig-1))*(2/Pi*lx/R+1.8642)*k;

C5(sig,lx,k)=8.008*C2(sig)*exp(-2*sqrt(lx/R))*lx/R*k;

C6(sig,lx)={local(k);k=k0(sig,lx);return(2.002*2^(5-2*sig)*C1(sig)+C3(sig,lx,k)+C4(sig,lx,k)+C5(sig,lx,k))};

e0(sig,lx)=C6(sig,lx)*exp((10-16*sig)/3*sqrt(lx/R))*(lx/R)^((5-2*sig)/2);



