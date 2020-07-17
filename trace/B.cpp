#include "stdio.h"
#include "stdlib.h"
#include "inttypes.h"
#include "../includes/int_double12.0.h"

#define N_PI_POWERS (3)
#define N_A_POWERS (6)
double a;

int_double pi_powers[N_PI_POWERS+1];
int_double a_powers[N_A_POWERS+1];
int_double b2,c,c_bit;

void setup()
{
  pi_powers[0]=1;
  a_powers[0]=1;
  for(uint64_t i=0;i<N_PI_POWERS;i++)
    pi_powers[i+1]=pi_powers[i]*d_pi;
  for(uint64_t i=0;i<N_A_POWERS;i++)
    a_powers[i+1]=a_powers[i]*a;
  b2=(6.0*d_pi*d_pi-1)/2.0;
  c=40.0*pi_powers[3]*a_powers[3]/(11*pi_powers[2]*a_powers[2]*b2-5.0);
  c_bit=c/(240.0*pi_powers[3]*a_powers[6]);
}

int_double beta_0_a(int_double *t_powers) // beta for t in [0,a]                                                                                                                                               
{
  return(2*c_bit*(pi_powers[2]*b2*(33.0*a_powers[5]-30.0*a_powers[3]*t_powers[2]+15.0*a*t_powers[4]-5*t_powers[5])-15.0*a_powers[3]+45.0*a*t_powers[2]-25.0*t_powers[3]));
}

int_double beta_a_2a(int_double *t_powers)
{
  return(c_bit*(pi_powers[2]*b2*(51.0*a_powers[5]+75.0*a_powers[4]*t_powers[1]-210.0*a_powers[3]*t_powers[2]+150.0*a_powers[2]*t_powers[3]-45.0*a*t_powers[4]+5.0*t_powers[5])-105.0*a_powers[3]+225.0*a_power\
s[2]*t_powers[1]-135.0*a*t_powers[2]+25.0*t_powers[3]));
}

int_double beta_2a_3a(int_double *t_powers)
{
  int_double t1=3*a-t_powers[1];
  t1=t1*t1*t1;
  return(c_bit*t1*(pi_powers[2]*b2*(9.0*a_powers[2]-6.0*a*t_powers[1]+t_powers[2])+5.0));
}

int_double eval_integral(uint64_t inv_step_size)
{
  double step_size=1.0/inv_step_size;
  // do from 0 to a                                                                                                                                                                                            
  double t=step_size;
  int_double t_powers[6];
  int_double res=0.0;
  for(uint64_t i=0;i<5;i++)
    t_powers[i+1]=t_powers[i]*t;

  int_double left=0.25;
  int_double right=(cosh(t*d_pi)-1)*beta_0_a(t_powers)/(2*(pi_powers[2]*t_powers[2]));
  int_double term;
  term.left=min(left.left,right.left);
  term.right=min(left.right,right.right);
  res+=term*step_size;
  while(t<a)
    {
      t+=step_size;
      left=right;
      for(uint64_t i=0;i<5;i++)
        t_powers[i+1]=t_powers[i]*t;
      right=(cosh(t*d_pi)-1)*beta_0_a(t_powers)/(2*(pi_powers[2]*t_powers[2]));
      term.left=min(left.left,right.left);
      term.right=min(left.right,right.right);
      res+=term*step_size;
    }
  while(t<2*a)
    {
      t+=step_size;
      left=right;
      for(uint64_t i=0;i<5;i++)
        t_powers[i+1]=t_powers[i]*t;
      right=(cosh(t*d_pi)-1)*beta_a_2a(t_powers)/(2*(pi_powers[2]*t_powers[2]));
      term.left=min(left.left,right.left);
      term.right=min(left.right,right.right);
      res+=term*step_size;
    }
  while(t<3*a)
    {
      t+=step_size;
      left=right;
      for(uint64_t i=0;i<5;i++)
        t_powers[i+1]=t_powers[i]*t;
      right=(cosh(t*d_pi)-1)*beta_2a_3a(t_powers)/(2*(pi_powers[2]*t_powers[2]));
      term.left=min(left.left,right.left);
      term.right=min(left.right,right.right);
      res+=term*step_size;
    }

  return(2*res);
}

int main(int argc, char **argv)
{
  if(argc!=3)
    {
      printf("Usage:- %s < a > < -log_2 (step_size) >\n",argv[0]);
      exit(0);
    }
  a=atof(argv[1]);
  uint64_t num_steps=1<<atol(argv[2]);
  _fpu_rndd();
  setup();
  printf("a=%10.8e\n",a);
  print_int_double_str("b^2=",b2);
  print_int_double_str("c=",c);
  int_double t_powers[6];
  t_powers[0]=1;
  for(uint64_t i=0;i<5;i++)
    t_powers[i+1]=t_powers[i]*0.5*a;
  print_int_double_str("beta(0.5a)=",beta_0_a(t_powers));
  for(uint64_t i=0;i<5;i++)
    t_powers[i+1]=t_powers[i]*1.5*a;
  print_int_double_str("beta(1.5a)=",beta_a_2a(t_powers));
  for(uint64_t i=0;i<5;i++)
    t_powers[i+1]=t_powers[i]*2.5*a;
  print_int_double_str("beta(2.5a)=",beta_2a_3a(t_powers));
  print_int_double_str("int -3a..3a =",eval_integral(num_steps));
  return(0);
}
