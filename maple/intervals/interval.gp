al=2;ah=3;bl=1;bh=2;
cl=3.5;ch=4;dl=0.5;dh=1.5;

rl2=(al^2+bl^2)*(cl^2+dl^2);
rh2=(ah^2+bh^2)*(ch^2+dh^2);

print("r^2=[",rl2,",",rh2,"]");
print("r=[",sqrt(rl2),",",sqrt(rh2),"]");

thetal=atan(bl/ah)+atan(dl/ch);
thetah=atan(bh/al)+atan(dh/cl);
print("theta=[",thetal,",",thetah,"]");
tl2=tan(thetal)^2;
th2=tan(thetah)^2;

xl2=(rh2*tl2-rl2)/(th2*tl2-1);
xh2=rh2-xl2*th2;
yl2=rl2-xl2;
yh2=rh2-xh2;

print("[",sqrt(xl2),",",sqrt(xh2),"]+[",sqrt(yl2),",",sqrt(yh2),"]i");

print(14*thetal);