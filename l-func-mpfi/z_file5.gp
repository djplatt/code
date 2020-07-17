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

zeta_file(101,512,14440,0.078125,"$WORKDIR/z_file_14440_14480.dat")
zeta_file(101,512,14480,0.078125,"$WORKDIR/z_file_14480_14520.dat")
zeta_file(101,512,14520,0.078125,"$WORKDIR/z_file_14520_14560.dat")
zeta_file(101,512,14560,0.078125,"$WORKDIR/z_file_14560_14600.dat")
zeta_file(101,512,14600,0.078125,"$WORKDIR/z_file_14600_14640.dat")
zeta_file(101,512,14640,0.078125,"$WORKDIR/z_file_14640_14680.dat")
zeta_file(101,512,14680,0.078125,"$WORKDIR/z_file_14680_14720.dat")
zeta_file(101,512,14720,0.078125,"$WORKDIR/z_file_14720_14760.dat")

quit
