/* produce G files for Odd Dirichlet Character L-function */

\p 200

eps=1e-100;

g(u)=2*exp(3*u/2-Pi*exp(2*u));

gn(u,n)=if(n==0,g(u),(g(u+eps)-g(u-eps))/(2*eps));

