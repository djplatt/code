/* generate coefficients for quadratic character mod 11 */
N=11;M=100;
print(N); /* conductor */
print(1); /* its self dual */
print(M); /* M */
for(m=1,M,print(kronecker(m,N)));
quit;
