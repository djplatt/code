x0=10^22;
x=x0-273747*2^32+1;

zeros_error=0.249;
sieve_error=0.25;

a=201467263477358935023.594; /* G(1)-G(1/2) */
a_err=0.001;
b=1.61227245535e8; /* n*G(inf)-sigma G(rho) - G(14) */
b_err=0.001;
c=-4.216110143713e9; /* G(14)-G(1/2) */
c_err=0.001
d=9.8043975189e7; /* G(inf)-G(14) */
d_err=0.001;
e=2059488903.355; /* Pi*(x)-Pi(x) */
e_err=0.001;
f=193818.795+1.439; /* sum phi(p)+sum phi(p^n) */
f_err=0.003+0.001;
log2=floor(log(2)*1000)/1000;
log2_err=0.001;

primes_to_10_22=23209771833641; /* from sieving process */

calc_pi_low=a+2*b-c-c_err-3*d-3*d_err-log2-log2_err-e-e_err+f-zeros_error-sieve_error+primes_to_10_22;
calc_pi_hi=a+a_err+2*b+2*b_err-c-3*d-log2-e+f+f_err+zeros_error+sieve_error+primes_to_10_22;

pi_x=201467286689315906290; /* from prime pages */

print("Calc Pi(",x0,") in [",calc_pi_low,",",calc_pi_hi,"]");
print(" Act Pi(",x0,") in [",pi_x,",",pi_x,"]");
print("Diff=[",calc_pi_low-pi_x,",",calc_pi_hi-pi_x,"]");
