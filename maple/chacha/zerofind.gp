zerofind(q,chi)={L=lfuninit(Mod(chi,q),[1/4,1/4,12],2);z=0;for(n=1,10,z=z-lfun(L,z,1)/lfun(L,z,2));if(abs(lfun(L,z,1))>1e-30,printf("q = %d chi = %d failed.\n",q,chi),printf("q = %d chi = %d zero near %f + %f i\n",q,chi,real(z),imag(z)))};

read("zerofind_test.gp");
quit;
