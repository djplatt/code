max_odd_q=7.5e5;
max_even_q=1.5e6;
max_odd_phi=150000;
max_even_phi=300000;
odd_t0_phi=3.75e7;
even_t0_phi=7.5e7;
t0_margin=50;
/*
print("Number of odd q to re-compute=",sum(q=3,300000,if(q%2==1,if(eulerphi(q)<max_odd_phi,if(1e8/q<(odd_t0_phi/eulerphi(q)+50),print(q);1)))));

print("List of odd q to re-compute with t0>=10000 ",for(q=3,300000,if(q%2==1,if(eulerphi(q)<max_odd_phi,t0=odd_t0_phi/eulerphi(q)+50;if((t0>1e8/q)&(t0>10000),print(q))))));

print("Cost of re-computing=",sum(q=3,300000,if(q%2==1,if(eulerphi(q)<max_odd_phi,if(1e8/q<(odd_t0_phi/eulerphi(q)+50),eulerphi(q)/1e8*(odd_t0_phi/eulerphi(q)+50))))));


print("Number of odd q to do from scratch=",sum(q=300000,max_odd_q,if(q%2==1,if(eulerphi(q)<max_odd_phi,print(q);1))));

print("Number of odd q to do from scratch=",sum(q=300000,max_odd_q,if(q%2==1,if(eulerphi(q)<max_odd_phi,print(odd_t0_phi/eulerphi(q)+50);1))));
*/
print("Number of even q to compute=",sum(q=4,1500000,if(q%4==0,ph=eulerphi(q);if(ph<max_even_phi,1))));

print("Number of even q with t0>=10000 ",sum(q=4,max_even_q,if(q%4==0,ph=eulerphi(q);if(ph<max_even_phi,t0=even_t0_phi/ph+50;if(t0>10000,1)))));

print("Cost even q with t0<10000 ",sum(q=4,max_even_q,if(q%4==0,ph=eulerphi(q);if(ph<max_even_phi,t0=even_t0_phi/ph+50;if(t0<1000,t0*ph/1e8)))));

