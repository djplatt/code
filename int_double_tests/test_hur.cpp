#include <stdio.h>
#include <stdlib.h>
#include "../includes/int_double11.0.h"

int main()
{
  int_complex hur;
  int worst_err=-20,r_err,i_err,m_err;
  _fpu_rndd();
  for(int i=1;i<1000;i++)
    for(int t=100;t<=10000;t+=100)
      {
	hur=hurwitz(int_complex(int_double(0.5),int_double(t)),int_double(i)/1000.0);
	r_err=abs_error(hur.real);
	i_err=abs_error(hur.imag);
	m_err=max(r_err,i_err);
	if(m_err>worst_err)
	  {
	    printf("for t=%d, alpha=%10.8e, absolute errors are %d and %d\n",t,((double)i)/1000.0,r_err,i_err);
	    print_int_complex_str("",hur);
	    worst_err=m_err;
	  }
      }
  return(0);
}
