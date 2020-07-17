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

zeta_file(101,512,15080,0.078125,"$WORKDIR/z_file_15080_15120.dat")
zeta_file(101,512,15120,0.078125,"$WORKDIR/z_file_15120_15160.dat")
zeta_file(101,512,15160,0.078125,"$WORKDIR/z_file_15160_15200.dat")
zeta_file(101,512,15200,0.078125,"$WORKDIR/z_file_15200_15240.dat")
zeta_file(101,512,15240,0.078125,"$WORKDIR/z_file_15240_15280.dat")
zeta_file(101,512,15280,0.078125,"$WORKDIR/z_file_15280_15320.dat")
zeta_file(101,512,15320,0.078125,"$WORKDIR/z_file_15320_15360.dat")
zeta_file(101,512,15360,0.078125,"$WORKDIR/z_file_15360_15400.dat")

quit
