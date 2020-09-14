print_int(str,a,b)=print(str," [",a,",",b,"]");

x=10^22;
x0=x-273747*2^32+1;
prec=0.001;

pi_x=201467286689315906290; /* primepages */


log2_left=floor(log(2)*1000)/1000.0;
log2_right=log2_left+prec;

G_1_half_left=201467263477358935023.594;
G_1_half_right=G_1_half_left+prec;

G_half_14i_left=-4216110143.713;
G_half_14i_right=G_half_14i_left+prec;

G_inf_14i_left=98043975.189;
G_inf_14i_right=G_inf_14i_left+prec;

sigma_G_left=161227245.535;
sigma_G_right=sigma_G_left+prec;

sigma_phi_left=193818.795; /* sigma p */
sigma_phi_right=193818.799;

sigma_phin_left=1.439; /* sigma p^n n>1 */
sigma_phin_right=sigma_phin_left+prec;

pi_star_less_pi_left=2059488903.355;
pi_star_less_pi_right=pi_star_less_pi_left+prec;

prime_count=23209771833641;

pi_x0_act=pi_x-prime_count;

pi_star_x0_act=pi_x0_act+pi_star_less_pi_left;

pi_star_x0_left=G_1_half_left+G_half_14i_right+3*G_inf_14i_left+2*sigma_G_left-log2_right+sigma_phi_left+sigma_phin_left;


ls=G_1_half_left;
rs=G_1_half_right;
print_int("We have G(1)-G(1/2) in",ls,rs);
l=G_half_14i_left;
r=G_half_14i_right;
print_int("We have G(1/2)-G(1/2+14i) in",l,r);
ls=ls+l;rs=rs+r;
print_int("sum so far in",ls,rs);
l=2*sigma_G_left;
r=2*sigma_G_right;
print_int("We have -2*sigma G-2*G(1/2+14*i) in",l,r);
ls=ls+l;rs=rs+r; 
print_int("sum so far in",ls,rs);
l=3*G_inf_14i_left;
r=3*G_inf_14i_right;
print_int("We have 3*G(1/2+14i) in",l,r);
ls=ls+l;rs=rs+r; 
print_int("sum so far in",ls,rs);
l=sigma_phi_left+sigma_phin_left;
r=sigma_phi_right+sigma_phin_right;
print_int("We have sigma 1/n [chi_x(p^n)-phi(p^n)] in",l,r);
ls=ls+l;rs=rs+r; 
print_int("sum so far in",ls,rs);
l=-log2_right;
r=-log2_left;
print_int("We have -log(2) in",l,r);
ls=ls+l;rs=rs+r; 
print_int("Pi*(x0) is in",ls,rs);
l=-pi_star_less_pi_right;
r=-pi_star_less_pi_left;
print_int("sum_{n>1} 1/n pi(x_0^{1/n}) in",l,r);
ls=ls+l;rs=rs+r;
print_int("Pi(x0) in ",ls,rs);
print("Pi(10^22)-Pi(x0)=",prime_count);
ls=ls+prime_count;
rs=rs+prime_count;
print_int("Pi(10^22) in",ls,rs);
print_int("error in",-0.5,0.5);
ls=ls-0.5;
rs=rs+0.5;
print_int("Pi(10^22) in",ls,rs);
print(    "Pi(10^22) =   ",pi_x);

