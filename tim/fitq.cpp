//
// fitq.cpp
//
// DJ Platt January 2014
//
// find a Q>=Q0 such that
// |zeta(sigma+it)|<a|it+sigma+Q|^1/6 log|it+sigma+Q|
// where the zeta values are in infile

#include "inttypes.h"
#include "../includes/int_double12.0.h"

int_double calcrhs(double q, int_double t2, int_double sigma, int_double a)
{
  int_double s=sqrt(t2+sqr(sigma+q));
  return(a*pow(s,int_double(1)/6)*log(s));
}

//
// find a q s.t. q>=q0 and z<=|q+1/2+It|^(1/6)*log(|q+1/2+it|)
//
// do at most this many iterations when binary chopping.
//
#define NUM_ITS (50)
double fitq(double q0, double t, double z, int_double sigma, int_double a)
{
  int_double t2=sqr(int_double(t));
  int_double rhs=calcrhs(q0,t2,sigma,a);
  if(z<=rhs.left) // q0 is already good enough
    return(q0);
  double newq=q0*1.01; // go up by 1% each time
  int_double new_rhs=calcrhs(newq,t2,sigma,a);
  while(z>new_rhs.left)
    {
      q0=newq;
      newq*=1.01;
      new_rhs=calcrhs(newq,t2,sigma,a);
    }
  // now newq is too big and q0 is too small
  // binary chop on [q0,newq]
  for(uint64_t i=0;i<NUM_ITS;i++)
    {
      double midq=(q0+newq)/2.0;
      if((midq<=q0)||(midq>=newq))
	return(newq);
      int_double mid_rhs=calcrhs(midq,t2,sigma,a);
      if(z>=mid_rhs.left)
	q0=midq;
      else
	newq=midq;
    }
  return(newq);
}

int main(int argc, char **argv)
{
  printf("Command Line:- %s ",argv[0]);
  for(uint64_t i=1;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=4)
    {
      printf("Usage:- %s <a*10000> <infile> <Q0>.\n",argv[0]);
      exit(0);
    }

  // set up the SSE system for rounding and initialise constants
  _fpu_rndd();

  int_double a=atol(argv[1]);
  a/=10000.0;
  print_int_double_str("Using a in ",a);
  FILE *infile=fopen(argv[2],"rb");
  if(!infile)
    {
      printf("Failed to open file %s for binary input. Exiting.\n",argv[2]);
      exit(0);
    }
  double q=atof(argv[3]);

  double t,z,sigma=0.5,worst_t=0.0,q1,worst_z=0.0;
  while(fread(&t,sizeof(double),1,infile)==1) // just read to end of file
    {
      if(fread(&z,sizeof(double),1,infile)!=1)
	{
	  printf("Error reading largest so far for t=%40.38e. Exiting.\n",t);
	  exit(0);
	}
      //printf("Doing t=%40.38e\n",t);
      q1=fitq(q,t,z,sigma,a);
      if(q1>q) // we have to increase q
	{
	  //printf("increasing q to %40.38e to cope with z=%40.38e at t=%40.38e.\n",q1,z,t);
	  worst_z=z;
	  worst_t=t;
	  q=q1;
	}
    }

  printf("Setting q to %40.38e will suffice at t=%40.38e where |z|=%40.38e.\n",q,worst_t,worst_z);
  return(0);
}
