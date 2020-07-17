#include <stdio.h>
#include "/cygdrive/c/home/gmp-4.2.2/gmp.h"
#include <iostream>
#define max(a,b) (a>=b?a:b)


class mpf_int
{
public:
  int prec;
  mpf_t left;
  mpf_t right;

  /*
  // the destructor
  inline ~mpf_int()
  {
    mpf_clear(left);
    mpf_clear(right);
  }
  */

  // constructor with default precision
  inline mpf_int()
  {
    prec=mpf_get_default_prec();
    mpf_init2(left,prec);
    mpf_init2(right,prec);
  };

  // constructor with precision supplied
  inline mpf_int(int this_prec)
  {
    prec=this_prec;
    mpf_init2(left,prec);
    mpf_init2(right,prec);
  };

  friend mpf_int operator + (const mpf_int &lhs, const mpf_int &rhs);

  friend mpf_int operator + (const mpf_int &lhs, const unsigned int rhs);

  friend mpf_int operator + (const mpf_int &lhs, const double rhs);


};

void print_mpf_int(const mpf_int &y)
{
  double l,r;
  l=mpf_get_d(y.left);
  r=-mpf_get_d(y.right);
  printf("[%e,%e]",l,r);
}


// add two mpf_int's
mpf_int operator + (const mpf_int &lhs, const mpf_int &rhs)
{
  int max_prec=max(lhs.prec,rhs.prec);
  mpf_int res(max_prec);
  mpf_add(res.left,lhs.left,rhs.left);
  mpf_add(res.right,lhs.right,rhs.right);
  return(res);
}

// add an mpf_int to 
mpf_int operator + (const mpf_int &lhs, const unsigned int rhs)
{
  mpf_int res(lhs.prec);
  mpf_add_ui(res.left,lhs.left,rhs);
  mpf_sub_ui(res.right,lhs.right,rhs);
  return(res);
}

mpf_int operator + (const mpf_int &lhs, const double rhs)
{
  printf("In + with ");
  print_mpf_int(lhs);
  printf(" and %e\n",rhs);
  mpf_int res(lhs.prec);
  mpf_set_d(res.left,rhs);
  mpf_set_d(res.right,-rhs);
  mpf_add(res.left,lhs.left,res.left);
  mpf_add(res.right,lhs.right,res.right);
  printf("Returning ");
  print_mpf_int(res);
  printf("\n");
  return(res);
}


int main()
{
  mpf_set_default_prec(100);
  mpf_int x,y(200); // x will have 100 ish bits, y 200 ish
  mpf_set_d(x.left,10.0);
  mpf_set_d(x.right,-11.0);
  y=x+(double)100.4;
  print_mpf_int(y);
  y=y+100.4;
  print_mpf_int(x);
  print_mpf_int(y);
  return(0);
}
