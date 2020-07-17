#include <stdio.h>
#include <stdlib.h>
#include "../includes/int_double12.0.h"
#define POS (0)
#define NEG (1)
#define TO_ZERO (2)
#define FROM_ZERO (3)
#define CROSS_ZERO (4)

#define my_printc(str,val) {printf(str);printf(": ");print_int_complex(val);printf("\n");}
#define my_print(str,val) {printf(str);printf(": ");print_int_double(val);printf("\n");}

double rand_double(char type)
{
	double x;
	switch(type)
	{
	case POS:
		x=rand()/RAND_MAX;
		return(x);
	case NEG:
		x=rand()/RAND_MAX;
		return(-x);
	case TO_ZERO:
		return(0.0);
	default:
		printf("Fatal error in rand_double. Exiting.\n");
		exit(0);
		return(x);
	}
}


int_double rand_int_double(char type)
{
	int_double x;

	switch(type)
	{
	case POS:
		x.left=rand();
		x.right=-x.left-1;
		x=x/RAND_MAX;
		return(x);
	case NEG:
		x.left=rand();
		x.right=-x.left-1;
		x=x/RAND_MAX;
		return(-x);
	case TO_ZERO:
		x.left=-rand();
		x.right=0.0;
		x=x/RAND_MAX;
		return(x);
	case FROM_ZERO:
	x.left=0.0;
		x.right=-rand();
		x=x/RAND_MAX;
		return(x);
	case CROSS_ZERO:
		x.left=-rand();
		x.right=-rand();
		x=x/RAND_MAX;
		return(x);
	default:
		printf("Fatal error in rand_int_double. Exiting.\n");
		exit(0);
		return(x);
	}
}

// ((x*y)/x)/y)=1
bool testa (unsigned int n)
{
	int_double x,y,z;
	unsigned int i;
	char type1,type2;

	for(type1=POS;type1<=NEG;type1++)
		for(type2=POS;type2<=NEG;type2++)
			for(i=0;i<n;i++)
			{
				x=rand_int_double(type1);
				y=rand_int_double(type2);
				z=x*y;
				z=z/x;
				z=z/y;
				z=z-1;
				if(!contains_zero(z))
				{
					my_print("x",x);
					my_print("y",y);
					my_print("x*y",z=x*y);
					my_print("x*y/x",z=z/x);
					my_print("x/y/x/y",z=z/y);
					return(false);
				}
			}
	return(true);
}



//(x+y)(x-y)=x^2-y^2
bool testb (unsigned int n)
{
	int_double x,y;
	unsigned int i;
	char type1,type2,type3,type4;

	for(type1=POS;type1<=CROSS_ZERO;type1++)
		for(type2=POS;type2<=CROSS_ZERO;type2++)
					for(i=0;i<n;i++)
					{
						x=rand_int_double(type1);
						y=rand_int_double(type2);
						if(!contains_zero(x*x-y*y-(x+y)*(x-y)))
						{
							printf("Test number %d.\n",i);
							my_print("x",x);
							my_print("y",y);
							my_print("x*x",x*x);
							my_print("y*y",y*y);
							my_print("x*x-y*y",x*x-y*y);
							my_print("x+y",x+y);
							my_print("x-y",x-y);
							my_print("(x+y)*(x-y)",(x+y)*(x-y));
							return(false);
						}
					}

	return(true);
}

//(x+y)(x-y)=x^2-y^2 with y a double
bool testc (unsigned int n)
{
	int_double x;
	double y;
	unsigned int i;
	char type1,type2,type3,type4;

	for(type1=POS;type1<=CROSS_ZERO;type1++)
		for(type2=POS;type2<=TO_ZERO;type2++)
					for(i=0;i<n;i++)
					{
						x=rand_int_double(type1);
						y=rand_double(type2);
						if(!contains_zero(x*x-y*y-(x+y)*(x-y)))
						{
							printf("Test number %d.\n",i);
							my_print("x",x);
							printf("y %20.18e\n",y);
							my_print("x*x",x*x);
							my_print("y*y",y*y);
							my_print("x*x-y*y",x*x-y*y);
							my_print("x+y",x+y);
							my_print("x-y",x-y);
							my_print("(x+y)*(x-y)",(x+y)*(x-y));
							return(false);
						}
					}

	return(true);
}


// (x^2-y^2)=(x+y)(x-y)
bool test1(unsigned int n)
{
	int_complex x,y;
	unsigned int i;
	char type1,type2,type3,type4;

	for(type1=POS;type1<=CROSS_ZERO;type1++)
		for(type2=POS;type2<=CROSS_ZERO;type2++)
			for(type3=POS;type3<=CROSS_ZERO;type3++)
				for(type4=POS;type4<=CROSS_ZERO;type4++)
					for(i=0;i<n;i++)
					{
						x=int_complex(rand_int_double(type1),rand_int_double(type2));
						y=int_complex(rand_int_double(type1),rand_int_double(type2));
						if(!contains_zero(x*x-y*y-(x+y)*(x-y)))
						{
							printf("Test number %d.\n",i);
							my_printc("x",x);
							my_printc("y",y);
							my_print("R(x*x)",x.real*x.real-x.imag*x.imag);
							my_print("I(x*x)",x.real*x.imag*2);
							my_printc("x*x",x*x);
							my_printc("y*y",y*y);
							my_printc("x*x-y*y",x*x-y*y);
							my_printc("x+y",x+y);
							my_printc("x-y",x-y);
							my_printc("(x+y)*(x-y)",(x+y)*(x-y));
							return(false);
						}
					}

	return(true);
}

// (x^2-y^2)/(x+y)=(x-y)
bool test2(unsigned int n)
{
	int_complex x,y,x_plus_y;
	unsigned int i;
	char type1,type2,type3,type4;
	for(type1=POS;type1<=CROSS_ZERO;type1++)
		for(type2=POS;type2<=CROSS_ZERO;type2++)
			for(type3=POS;type3<=CROSS_ZERO;type3++)
				for(type4=POS;type4<=CROSS_ZERO;type4++)
					for(i=0;i<n;i++)
					{
						x=int_complex(rand_int_double(type1),rand_int_double(type2));
						y=int_complex(rand_int_double(type1),rand_int_double(type2));
						x_plus_y=x+y;
						if(!contains_zero(x_plus_y))
							if(!contains_zero((x*x-y*y)/x_plus_y-(x-y)))
								return(false);
					}

	return(true);
}

//sin^2(x)+cos^2(x)=1
bool test3(unsigned int n)
{
	int_double x,cos_x,sin_x;
	unsigned int i;
	char type1;
	for(x.left=0.0;x.left<M_PI/2;x.left+=1.0/(double) n)
	  {
	    x.right=-x.left-0.00000001;
	    sin_cos(x,&sin_x,&cos_x);
	    if(!contains_zero(sqr(sin_x)+sqr(cos_x)-1))
	      {
		printf("x: ");
		print_int_double(x);
		printf("\nsin(x): ");
		print_int_double(sin_x);
		printf("\ncos(x): ");
		print_int_double(cos_x);
		printf("\n");
		return(false);
	      }
	  }
	return(true);
}

//log(exp(x))=x
bool test4 (unsigned int n)
{
	int_double x;
	unsigned int i;
	char type1;
	for(type1=POS;type1<=CROSS_ZERO;type1++)
		for(i=0;i<n;i++)
		{
			x=rand_int_double(type1);
			if(!contains_zero(log(exp(x))-x))
				return(false);
		}
		return(true);
}
	
//exp(log(x))=x
bool test5 (unsigned int n)
{
    int_double x,lnx,explnx;
	unsigned int i;
	for(i=0;i<n;i++)
	{
		x=rand_int_double(POS);
		lnx=log(x);
		explnx=exp(log(x));
		if(!contains_zero(explnx-x))
		{
		    printf("i: %d\nx: ",i);
		    print_int_double(x);
		    printf("\nlog(x): ");
		    print_int_double(log(x));
		    printf("\nlog(-x.right): %20.18e\nnextafter(log(-x.right)): %20.18e",
			   log(-x.right),nextafter(log(-x.right)));
		    printf("\nexp(log(x)): ");
		    print_int_double(exp(log(x)));
		    printf("\nexp(lnx.left): %20.18e\n",exp(lnx.left));
		    printf("exp(-lnx.right): %20.18e\n",exp(-lnx.right));

		    printf("\nexp(log(x))-x: ");
		    print_int_double(exp(log(x))-x);
      		    return(false);
		}
	}
	return(true);
}

//atan2(sin(x),cos(x))=x
bool test6 (unsigned int n)
{
  double theta0;
  int_double theta,si,co,atn;
  for(theta0=-3.14/2.0;theta0<=3.13;theta0+=1.0/(double) n)
    {
      theta=int_double(theta0,theta0+1.0/100.0);
      sin_cos(theta,&si,&co);
      atn=atan2(si,co);
      if(!(contains_zero(atn-theta)))
	{
	  print_int_double_str("theta=",theta);
	  print_int_double_str("sin  =",si);
	  print_int_double_str("cos  =",co);
	  print_int_double_str("atan =",atn);
	  return(false);
	}
    }
  return(true);
}

// cos(2x)=2cos^2(x)-1
bool test7 (unsigned int n)
{
  double theta0;
  int_double theta,si,co1,co2;
  for(theta0=-3.14;theta0<=3.13;theta0+=1.0/(double) n)
    {
      theta=int_double(theta0,theta0+1.0/100.0);
      sin_cos(theta*2,&si,&co1);
      sin_cos(theta,&si,&co2);
      if(!(contains_zero(co1+1.0-2*sqr(co2))))
	{
	  print_int_double_str("theta=",theta);
	  print_int_double_str("co1  =",co1);
	  print_int_double_str("co2  =",co2);
	  return(false);
	}
    }
  return(true);
}

bool test8 (unsigned int n)
{
  int_double x;
  for(int i=0;i<n;i++)
    {
      x=rand_int_double(POS)*i;
      if(!contains_zero(sqr(sqrt(x))-x))
	return(false);
      x=rand_int_double(FROM_ZERO)*i;
      if(!contains_zero(sqr(sqrt(x))-x))
	return(false);
    }
  return(true);
}

int_double intersection(const int_double &x, const int_double &y)
{
  int_double res;
  //print_int_double_str("x",x);
  //print_int_double_str("y",y);
  res.left=max(x.left,y.left);
  res.right=max(x.right,y.right);
  return(res);
}

int_complex fixup(const int_complex &z)
{
  int_double x;
  int_complex res;

  x=sqrt(d_one-sqr(z.real));
  if(z.imag.left<0)
    res.imag=-intersection(-z.imag,x);
  else
    res.imag=intersection(z.imag,x);

  x=sqrt(d_one-sqr(z.imag));
  if(z.real.left<0)
    res.real=-intersection(-z.real,x);
  else
    res.real=intersection(z.real,x);

  return(res);
}

int main()
{
	_fpu_rndd();
	printf("Testing....\n");
	print_int_double_str("Pi=",d_pi);
	print_int_double_str("Pi/2=",d_pi_2)
	print_int_double_str("2 Pi=",d_two_pi);
	print_int_double_str("sqrt(2)=",sqrt(int_double(2)));
	if(!testa(10000000))
	    printf("testa failed.\n");
	else
	    printf("testa passed.\n");

	if(!testb(10000000))
	    printf("testb failed.\n");
	else
	    printf("testb passed.\n");

	if(!testc(10000000))
	    printf("testc failed.\n");
	else
	    printf("testc passed.\n");
	return(0);
	if(!test1(100000))
	    printf("test1 failed.\n");
	else
	    printf("test1 passed.\n");
	if(!test2(100000))
	    printf("test2 failed.\n");
	else
	    printf("test2 passed.\n");
	if(!test3(100000))
	    printf("test3 failed.\n");
	else
	    printf("test3 passed.\n");
	if(!test4(100000))
	    printf("test4 failed.\n");
	else
	    printf("test4 passed.\n");
	if(!test5(100000))
	    printf("test5 failed.\n");
	else
	    printf("test5 passed.\n");
	
	if(!test6(100000))
	    printf("test6 failed.\n");
	else
	    printf("test6 passed.\n");
	if(!test7(100000))
	    printf("test7 failed.\n");
	else
	    printf("test7 passed.\n");
	if(!test8(100000))
	    printf("test8 failed.\n");
	else
	    printf("test8 passed.\n");

	return(0);
}

