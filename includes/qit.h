#ifndef QIT
#define QIT
#include "../includes/qit_struct.h"
#include "../includes/upsamdefs.h"

qit_t qit;

bool read_qit(FILE *qit_file, long unsigned int q)
{
  long int q_read;
  if(fread(&q_read,sizeof(long int),1,qit_file)!=1)
    {
      printf("Error reading qit file. Exiting.\n");
      exit(1);
    }
  if(q_read!=q)
    {
      printf("q in qit file (%ld) does not match required q (%lu). Exiting.\n",q_read,q);
      exit(1);
    }
  if(fread(&qit,sizeof(qit_t),1,qit_file)!=1)
    {
      printf("Error reading qit file. Exiting.\n");
      exit(1);
    }    
}

inline int_complex calc_qit(unsigned int n_dt)
{
  int_complex res=c_one;
  //printf("n_dt=%d\n",n_dt);
  for(int bit=0;bit<NUM_QIT_BITS;bit++,n_dt>>=1)
    if(n_dt&1)
      {
	res*=qit.bits[bit];
	//printf("bit %d is set\n",bit);
	//print_int_complex_str("res is now",res);
      }
  return(res);
}

inline int_complex q_it(unsigned int q, double t)
{
  double dn=t/(one_over_two_B/2.0);
  int n=dn;
  if((dn-(double) n)!=0.0)
    {
      printf("q_it called with non-integral multiple of 2B. Exiting.\n");
      exit(1);
    }
  return(calc_qit(n));
}
#endif
