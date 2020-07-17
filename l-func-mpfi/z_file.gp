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

zeta_file(101,512,8040,0.078125,"$HOME/data/oldlattice/z_file_8040_8080.dat")
zeta_file(101,512,8080,0.078125,"$HOME/data/oldlattice/z_file_8080_8120.dat")
zeta_file(101,512,8120,0.078125,"$HOME/data/oldlattice/z_file_8120_8160.dat")
zeta_file(101,512,8160,0.078125,"$HOME/data/oldlattice/z_file_8160_8200.dat")
zeta_file(101,512,8200,0.078125,"$HOME/data/oldlattice/z_file_8200_8240.dat")
zeta_file(101,512,8240,0.078125,"$HOME/data/oldlattice/z_file_8240_8280.dat")
zeta_file(101,512,8280,0.078125,"$HOME/data/oldlattice/z_file_8280_8320.dat")
zeta_file(101,512,8320,0.078125,"$HOME/data/oldlattice/z_file_8320_8360.dat")

quit
