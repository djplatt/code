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

zeta_file(101,512,15400,0.078125,"$WORKDIR/z_file_15400_15440.dat")
zeta_file(101,512,15440,0.078125,"$WORKDIR/z_file_15440_15480.dat")
zeta_file(101,512,15480,0.078125,"$WORKDIR/z_file_15480_15520.dat")
zeta_file(101,512,15520,0.078125,"$WORKDIR/z_file_15520_15560.dat")
zeta_file(101,512,15560,0.078125,"$WORKDIR/z_file_15560_15600.dat")
zeta_file(101,512,15600,0.078125,"$WORKDIR/z_file_15600_15640.dat")
zeta_file(101,512,15640,0.078125,"$WORKDIR/z_file_15640_15680.dat")
zeta_file(101,512,15680,0.078125,"$WORKDIR/z_file_15680_15720.dat")

quit
