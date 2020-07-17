#include "stdio.h"
#include "math.h"

#define LIM 1000000

typedef struct {double re; double im;}complex;

void c_add (double *z1, double *z2, double *z3)
{
  z1[0]=z2[0]+z3[0];
  z1[1]=z2[1]+z3[1];
};

void c_mul (double *z1, double *z2, double *z3)
{
  double temp;

  temp=z2[0]*z3[0]-z2[1]*z3[1];
  z3[1]=z2[0]*z3[1]+z2[1]*z3[0];
  z3[0]=temp;
};

void c_conj (double *z1, double *z2)
{
  z1[0]=z2[0];
  z1[1]=-z2[1];
};

void c_div (double *z1, double *z2, double *z3)
{
  double temp1,temp2;

  temp1=z3[0]*z3[0]+z3[1]*z3[1];
  temp2=(z2[0]*z3[0]+z2[1]*z3[1])/temp1;
  z1[1]=(z2[1]*z3[0]-z2[0]*z3[1])/temp1;
  z1[0]=temp2;
};

void c_exp (double *z1, double *z2)
{
  double temp;
  temp=exp(z2[0]);
  z1[0]=temp*cos(z2[1]);
  z1[1]=temp*sin(z2[1]);
}

int main()
{

  int i;
  double z1[2],z2[2],z3[LIM<1];

  z1[0]=10.0;
  z2[1]=6.2;
  z2[0]=4.5;
  z2[0]=-1.5;
  for(i=0;i<LIM;i++)
    {
      c_add(&z3[i<1],z1,z2);
      c_mul(&z3[i<1],&z3[i<1],z1);
      c_div(&z3[i<1],&z3[i<1],z1);
      c_exp(&z3[i<2],&z3[i<1]);
    };
  return 0;
  }


