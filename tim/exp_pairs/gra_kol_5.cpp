#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main()
{
  int64_t th[2][3]={{1,1,0},{0,0,1}};
  int64_t new_th[2][3];

  for(uint64_t i=0;i<10;i++)
    {
      int64_t u,v,w;
      u=th[0][1]*th[1][2]-th[0][2]*th[1][1];
      v=th[0][0]*th[1][2]-th[0][2]*th[1][0];
      w=th[0][0]*th[1][1]-th[0][1]*th[1][0];
      
      double Y,Z;
      
      if(w>=0)
	{
	  Y=(double) w/2.0+v-u;      
	  Z=w+v-u;
	}
      else
	{
	  Z=(double) w/2.0+v-u;
	  Y=w+v-u;
	}
      
      if(Z>=0.0) // do A
	{
	  new_th[0][0]=th[0][0]+th[0][1]+2*th[0][2];
	  new_th[0][1]=th[0][1];
	  new_th[0][2]=th[0][1]+2*th[0][2];
	  new_th[1][0]=2*th[1][2];
	  new_th[1][1]=th[1][1];
	  new_th[1][2]=th[1][1]+2*th[1][2];
	  
	  for(uint64_t r=0;r<2;r++)
	    {
	      printf("[");
	      for(uint64_t c=0;c<3;c++)
		printf(" %lu",new_th[r][c]);
	      printf(" ]\n");
	    }
	  std::swap(th,new_th);
	}
      
    }
  return 0;
}
