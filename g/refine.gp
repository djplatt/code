read("zeros1.gp");
refine(t)=
{
  local(d);
  d=-I/zeta'(1/2+I*t);
  t-=real(d*zeta(1/2+I*t));
  t-=real(d*zeta(1/2+I*t));
  t
}

go()=
{
  for(i=1,length(z),z[i]=refine(z[i]);if(i%1000==0,print(i)));
  write("zeros2.gp","z="z);
}
