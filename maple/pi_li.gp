/* returns a (low) estimate of the maximum number
   of primes in a window of width w
*/
max_primes(w)=w-(floor(w/2)+floor(w/3)-floor(w/6)+floor(w/5)-floor(w/10)-floor(w/15)+floor(w/30)-3);

max_primes1(w)=w*55296/323323;

A=9.682e147;
om_et=727.9513462695;
/*
y=solve(y=0,1e160,-max_primes1(y)+y/om_et-A);
print(y);
y1=solve(y=0,1e160,A-y/om_et);
print(y1);
*/