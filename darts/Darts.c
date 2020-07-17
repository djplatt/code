//
// Darts.c
// Created: 8th November 2009
//
// Version 1.0
#include <stdio.h>
#include "gmp.h"
#include "mpfr.h"
#define MIN(a, b) (((a)<(b)) ? (a) : (b))
#define mpz_inc(n) mpz_add_ui(n,n,1)
#define MAX_OUT (501) // playing 501 double out

mpz_t oneDarts[61]; // how many ways of scoring <n> with one dart
mpz_t allDarts[MAX_OUT+1]; // how many ways of finishing from <n>, last dart a double
mpz_t ways,tmp;
mpfr_t mpfr_tmp;

void InitOneDarts()
{
	int i;
	mpfr_set_default_prec(53);
	mpfr_init(mpfr_tmp);

	for(i=2; i<=MAX_OUT;i++)
		mpz_init(allDarts[i]);

	for(i=1; i<=60; ++i)
		mpz_init(oneDarts[i]);

	mpz_init(ways);
	mpz_init(tmp);

	for(i=1; i<=20; ++i)
	{
		mpz_inc(oneDarts[i]);     // single i
		mpz_inc(oneDarts[i+i]);   // double i
		mpz_inc(oneDarts[i+i+i]); // treble i
	}

	mpz_inc(oneDarts[25]); // bull
	mpz_inc(oneDarts[50]); // double bull
}

// assumes oneDarts has been filled
// and allDarts filled to left-1
void CountWays(int left)
{
	int i,j,dart;

	mpz_set_ui(ways,0); // count how many ways
	if(((left&1)==0)&&(left<=40)) // its a double
		mpz_inc(ways);            // so that's one way
	if(left==50)              // double bull
		mpz_inc(ways);        // so that's one way
	for(dart=1; dart<=MIN(left-2,60); ++dart) // score dart with first go
		                                      // must leave yourself at least 2
	{
		j = left - dart; // j needed to finish after hitting dart
		mpz_mul(tmp,oneDarts[dart],allDarts[j]);  // number of ways of hitting dart * number of ways out from j
		mpz_add(ways,ways,tmp);                   //
	}
	mpz_set(allDarts[left],ways);
}

#define FALSE (1==0)
#define TRUE (0==0)

// this answers Mike Harvey's question
void ThreeDarts()
{
	int gettable[181],i,j,k; // gettable[i] <=> can get i in three darts

	for(i=0;i<=180;i++)
		gettable[i]=FALSE; // don't know how to get anything yet

	for(i=0;i<=60;i++)
		if((i<=20)||(i==25)||(i==50)||(((i&1)==0)&&(i<=40))||(((i%3)==0)&&(i<=60))) // legal one dart scores
			for(j=0;j<=60;j++)
				if((j<=20)||(j==25)||(j==50)||(((j&1)==0)&&(j<=40))||(((j%3)==0)&&(j<=60))) // 2nd dart
					for(k=0;k<=60;k++)
						if((k<=20)||(k==25)||(k==50)||(((k&1)==0)&&(k<=40))||(((k%3)==0)&&(k<=60))) // 3rd dart
							gettable[i+j+k]=TRUE; // can get this

	for(i=0;gettable[i];i++); // skip to firt non-gettable score
	printf("%d is not gettable in three darts.\n",i); // print it
}




int main(unsigned long long int argc, char *argv[])
{
	int i;

	//ThreeDarts();
	//exit(0);
	InitOneDarts();
	mpz_set_ui(allDarts[2],1); // one way from 2, i.e. double 1

	for(i=3; i<=MAX_OUT; ++i)
		CountWays(i);
	for(i=2;i<=20;i++) // print the first 19
	{
		printf("Ways of getting out from %d (finishing on a double)=\n",i);
		mpz_out_str(NULL,10,allDarts[i]);printf("\n");
	}
	// print the 501 answer
	printf("Ways of getting out from %d (finishing on a double)=\n",MAX_OUT);
	mpz_out_str(NULL,10,allDarts[MAX_OUT]);printf("\n");
	return(0);
}

