#include <stdio.h>
#include <stdlib.h>
#include "../includes/int_double13.0.h"

int main()
{
	_fpu_rndd();
	/*
	int_double x,y,z;
	x=1.0;
	y=1e-30;

	z=x+y;
	print_int_double(z);printf("\n");

	z=x-y;
	print_int_double(z);printf("\n");

	z=-z;
	print_int_double(z);printf("\n");

	z=x*y;
	print_int_double(z);printf("\n");

	print_int_double(d_zero);printf("\n");
	print_int_double(d_neg_zero);printf("\n");
	print_int_double(d_neg_neg_zero);printf("\n");
	print_int_double(d_zero_zero);printf("\n");
	*/
	int_double bar(1e10),foo(1.1,1.2),bletch,bb;

	bletch=bar+100.0;
	//bletch=99.0+foo+bletch;
	bb=bar*bar;

	bar=bar+int_double(1e-100);

	printf("bar    %30.28f,%30.28f\n",left(bar),-right(bar));
	/*
	printf("foo    %30.28f,%30.28f\n",foo.left(),foo.right());
	printf("bletch %30.28f,%30.28f\n",bletch.left(),bletch.right());
	printf("bb     %30.28f,%30.28f\n",bb.left(),bb.right());
	bb=int_double(1)/int_double(3);
	printf("1/3    %30.28f,%30.28f\n",bb.left(),bb.right());
	bb=int_double(1)/int_double(-3);
	printf("-1/3    %30.28f,%30.28f\n",bb.left(),bb.right());
	bb=int_double(1)/3;
	printf("1/3    %30.28f,%30.28f\n",bb.left(),bb.right());
	bb=int_double(1)/(-3);
	printf("1/3    %30.28f,%30.28f\n",bb.left(),bb.right());
	*/
	if(int_double(1.1,1.2)>int_double(1.0,1.05))
	  printf("Good\n");
	if(int_double(1.1,1.2)>int_double(1.0,1.15))
	  printf("Bad.\n");
	int_double sinx,cosx;
	sin_cos(int_double(1,2),&sinx,&cosx);
	print_int_double_str("sin(1,2)=",sinx);
	printf("[%30.28e,%30.28e]\n",nextbefore(1e10),nextafter(1e10));
	return(0);
}

