\p 200
digits=180
ln10=log(10)

output(fname,x)=
{
   local(delta);
   if(x==0.0,delta=0,if(x<0.0,delta=10.0^(-digits+ceil(log(-x)/ln10)),\
                              delta=10.0^(-digits+ceil(log(x)/ln10))));
   write(fname,"[",x-delta,",",x+delta,"]");
}

zeta_file(N,num_s,im_start,step,fname)=
{
   local(z,z1,im,sm);

   write(fname,num_s);
/* N is max number of Taylor terms to use */
   write(fname,N);
   im=im_start;
   default(format,"e0.500");

   for(i=0,num_s-1,z1=0.5+im*I;\
             write(fname,im);\
             z=gamma(z1/2.0);\
             output(fname,real(z));\
             output(fname,imag(z));\
             z=gamma((z1+1.0)/2.0);\
             output(fname,real(z));\
             output(fname,imag(z));\
             forstep(re=0,N-1,1,\
                     z=zeta(z1);\
                     output(fname,real(z));\
                     output(fname,imag(z));\
                     z1=z1+1.0);\
            im=im+step);
}

zeta_file(101,512,14120,0.078125,"$WORKDIR/z_file_14120_14160.dat")
zeta_file(101,512,14160,0.078125,"$WORKDIR/z_file_14160_14200.dat")
zeta_file(101,512,14200,0.078125,"$WORKDIR/z_file_14200_14240.dat")
zeta_file(101,512,14240,0.078125,"$WORKDIR/z_file_14240_14280.dat")
zeta_file(101,512,14280,0.078125,"$WORKDIR/z_file_14280_14320.dat")
zeta_file(101,512,14320,0.078125,"$WORKDIR/z_file_14320_14360.dat")
zeta_file(101,512,14360,0.078125,"$WORKDIR/z_file_14360_14400.dat")
zeta_file(101,512,14400,0.078125,"$WORKDIR/z_file_14400_14440.dat")

quit
