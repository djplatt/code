#include "stdio.h"
#include "math.h"
#define N (30)
#define N_PRIMES (6)
#define MAX_PRIME (37)
int primes[N_PRIMES]={3,7,13,19,31,37};
int P,P2;
bool recips[MAX_PRIME*MAX_PRIME];

void set_recips()
{
	for(int i=0;i<P2;i++)
		recips[i]=false;
	for(int i=0;i<P2;i++)
	{
		long long int z=i;
		z=z*z*z;
		z=z%P2;
		recips[z]=true;
	}
}

bool check(int x, int y, int n)
{
	int z=x-y+n;
	z=z%P2;
	if(z<0)
		z+=P2;
	return(recips[z]);
}

int main()
{
	double ratio=1.0;
	for (int i=0;i<N_PRIMES;i++)
	{
		P=primes[i];
		P2=P*P;
		set_recips();

		int counter=0;
		for(int x=0;x<P2;x++)
			for(int y=0;y<P2;y++)
			{
				if(!check(x,y,N))
					counter++;
				if(!check(x,y,-N))
					counter++;}
			ratio*=1-((double) counter / (double) (2*P2*P2));
			printf("for p=%d a total of %d out of %d don't work n (%f).\n",P,counter,2*P2*P2,(double) counter/ (double)(2*P2*P2));
			printf("ratio now %e.\n",ratio);
	}
}
