#include "int_double4.0.cpp"
#include <iostream>
using namespace std;
int main()
{
  double t1,t2,t3,t4;
  int_double i1,rexp;

  _fpu_rndd;

  while(true)
    {

      cin >> t1 >> t2;

      if(t1==123.0)
	break;

      i1.left=t1;
      i1.right=-t2;
      printf("%20.18e %20.18e\n",i1.left,i1.right);
      printf("exp = ");
      print_int_double(exp(i1));
      printf("\nsqr = ");
      print_int_double(sqr(i1));
      printf("\nlog = ");
      print_int_double(log(i1));
      printf("\n\n");



    };
  return(0);
};

