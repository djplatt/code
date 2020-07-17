#include "../includes/int_double12.0.h"

int_double best_c;

void print_vars(int_double a,int_double b,int_double c, int_double d)
{
  print_int_double_str("a=",a);
  print_int_double_str("b=",b);
  print_int_double_str("c=",c);
  print_int_double_str("d=",d);
}

int_double intersect(int_double x, int_double y)
{
  int_double res;
  if(x.left>y.left)
    res.left=x.left;
  else
    res.left=y.left;
  if(x.right>y.right)
    res.right=x.right;
  else
    res.right=y.right;
  if(-res.right<res.left)
    {
      printf("Intersection of ");
      print_int_double(x);
      print_int_double_str(" has no intersection with ",y);
      exit(0);
    }
  return(res);
}

int_double c_from_ab(int_double a, int_double b)
{
  int_double res=(4*a*b+2)*(a+b-2*sqrt(a*b+1))+2*(a+b);
  if(res.left<4*b.left*a.left)
    res.left=4*b.left*a.left;
  return(res);
}

int_double dp(int_double a, int_double b)
{
  best_c=intersect(best_c,c_from_ab(a,b));
  return(a+b+best_c+2*a*b*best_c+2*sqrt((a*b+1)*(a*best_c+1)*(b*best_c+1)));
}

int_double five_thirds;

int_double d_from_c(int_double c)
{
  int_double low_d=pow(c,five_thirds);
  int_double high_d=sqr(c);
  int_double new_d;
  new_d.left=low_d.left;
  new_d.right=high_d.right;
  return(new_d);
}

int_double c_from_d(int_double d)
{
  int_double high_c=pow(d,0.6);
  int_double low_c=sqrt(d);
  return(int_double(low_c.left,-high_c.right));
}

#define r_min ((double) 1e7) // computation
#define AB_min r_min*r_min 
#define b_max (double) 4.3300001e23 // EFF 2013

int main()
{
  _fpu_rndd();
  int_double AB(AB_min,b_max*b_max);
  five_thirds=int_double(5);
  five_thirds/=3;
  int_double a(r_min/2,b_max);
  int_double b(r_min/2,b_max);
  int_double c(4*AB_min,2.03e70);
  int_double d(4*AB_min*c.left,-(a.right+b.right+c.right+2*a.right*b.right*c.right));
  d.right-=2*sqrt((a.right*b.right+1)*(a.right*c.right+1)*(b.right*c.right+1));
  print_vars(a,b,c,d);


  return(0);
}
