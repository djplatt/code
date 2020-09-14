al=2;ah=3;bl=1;bh=2;
cl=3.5;ch=4;dl=0.5;dh=1.5;
del=0.05;
min_del=del/10;
forstep(x=al,ah+min_del,del,forstep(y=bl,bh+min_del,del,print(x," ",y)));
forstep(x=cl,ch+min_del,del,forstep(y=dl,dh+min_del,del,print(x," ",y)));
forstep(x=al,ah+min_del,del,forstep(y=bl,bh+min_del,del,forstep(x1=cl,ch+min_del,del,forstep(y1=dl,dh+min_del,del,print(x*x1-y*y1," ",x*y1+y*x1)))));
quit

