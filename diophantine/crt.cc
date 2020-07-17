#include "stdio.h"
#include "malloc.h"
#include "math.h"
#define N (30)
#define N_PRIMES (6)
#define MAX_PRIME (37)
int Ps[N_PRIMES]={3,7,13,19,31,37};
int P2s[N_PRIMES]={9,49,169,19*19,31*31,37*37};

bool *recips[N_PRIMES];

void set_recips()
{
	for(int np=0;np<N_PRIMES;np++)
	{
		int P=Ps[np];
		int P2=P2s[np];
		recips[np]=(bool *) malloc(sizeof(bool)*P2);
		for(int i=0;i<P2;i++)	
			recips[np][i]=false;
		for(int i=0;i<P2;i++)
		{
			long long int z=i;
			z=z*z*z;
			z=z%P2;
			recips[np][z]=true;
		}
	}
}

bool check(int x, int y, int n, int P2)
{
	int z=x-y+n;
	z=z%P2;
	if(z<0)
		z+=P2;
	return(recips[z]);
}

// return solution to x=x1(m1),x=x2(m2) in [0,m1*m2)
inline int crt (int x1, int m1, int x2, int m2)
{
	return(0);
}

int main()
{
	set_recips();
}
