double QT;

inline int gcd (unsigned long int a, unsigned long int b)
// Euclid algorithm gcd
{
	unsigned int c;
	while(a!=0)
	{
		c=a;
		a=b%a;
		b=c;
	};
	return(b);
};

inline int co_prime(unsigned long int a, unsigned long int b)
{
	return(gcd(a,b)==1);
};


void set_QT_even(long unsigned int q)
{
  long unsigned int phi_q=1,d;
  for(d=2;d<q;d++)
    if(co_prime(q,d)) phi_q++;
  //printf("phi(%lu)=%lu.\n",q,phi_q);
  if(q&1)
    QT=(double) ceil(3.75e7/phi_q+50);
  else
    QT=(double) ceil(7.5e7/phi_q+50); // to eliminate Harald's middle arcs
  printf("t0 set to %f\n",QT);
  QT*=q;
  //printf("QT set to %f.\n",QT);
}
