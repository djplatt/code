/* compute valid intervals for integration of V
   so that pole of psi(1+2ir) at r=i/2 is outside
   |z|<=2
*/

comp_b(a,R)=solve(b=0.001,max(10,10*a),a-b+1/R/sin(atan(2/(a+b))));

max_B=10000;

printf("{0");
count=1;
b=comp_b(0,10.0); /* use |R|<=10.0 to be safe */

while(b<max_B,bb=floor(b*1024);printf(",%d",bb);b=comp_b(b,10.0);count++);

printf("}\n");

print(count);

quit;


