/* reads in the output from the LMFDB script and outputs it in an easier to handle format */
allocatemem(2^30);
x=readvec("temp.txt");
print(x[1]); /* conductor */
print(length(x[2])); /* number of gamma_r mus */
for(n=1,length(x[2]),print(x[2][n]));
if(length(x[3])!=0,print("This rep has Gamma_C's.");quit);
print(x[4]);
print(length(x[5]));
for(n=1,length(x[5]),for(m=1,length(x[5][n]),printf("%d ",x[5][n][m]));printf("-1\n"));
print(length(x)-5);
p=2;
for(n=6,length(x),print(p," ",x[n]);p=nextprime(p+1));
quit;
