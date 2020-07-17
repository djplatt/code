#include <iostream>
#include <float.h>

using namespace std;

unsigned long __nze[2]={0,0x80000000}; // -0.0

unsigned long __delta[2]={1,0};

double _nze = *(double*)__nze;

double _delta = *(double*)__delta;

//union ieee {double x; unsigned int i[2];};

inline double nextafter (double x)
{
  unsigned int *i=(unsigned int *) &x;




  if(i[1]&0x80000000)   // -ve
    {
      if(x==_nze)               // -0.0
	return(_delta);

      i[0]--;
      if(i[0]==0xffffffff)
	i[1]--;
      return(x);
    };

  if((i[1]&0x7ff00000)==0x7ff00000) // nan or +/-inf
    return(x);


  i[0]++;
  if(i[0]==0)
    i[1]++;
  return(x);
};

int main()
{

  double x;

  while(true)
    {
      cin >> x;
      printf("%20.18e\n",nextafter(x));
      if(x==123)
	exit(0);
    };

  return(0);
};
