#include <stdio.h>
#include <stdlib.h>
#include "../includes/int_double12.0.h"

int main()
{
	_fpu_rndd();

	int_double x,y;
	x=1.0;
	y=1e-30;
	x=x+y;
	print_int_double(x);printf("\n");
	x=x-y;
	print_int_double(x);printf("\n");

	printf("Testing....\n");
	return(0);
}

