print_int(str,a,b)=print(str," [",a,",",b,"]");

x=10^22;
x0=x-273747*2^32+1;
prec=0.001;

pi_x=201467286689315906290; /* primepages */


log2l=floor(log(2)*1000)/1000.0;
log2r=log2l+prec;

al=201467263477358935023.594;
ar=al+prec;

bl=-4216110143.713;
br=bl+prec;

cl=98043975.189;
cr=cl+prec;

dl=161227245.535;
dr=dl+prec;

el=193818.795; /* sigma p */
er=193818.799;

fl=1.439; /* sigma p^n n>1 */
fr=fl+prec;

gl=2059488903.355;
gr=gl+prec;

h=23209771833641;

errl=-0.5;
errr=0.5;

print("Lower bound = ",al-br+2*dl-3*cr+el+fl-gr+h+errl-log2r);
print("Upper bound = ",ar-bl+2*dr-3*cl+er+fr-gl+h+errr-log2l);