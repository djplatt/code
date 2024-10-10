#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int64_t gcd (int64_t a1, int64_t b1)
/* Euclid algorithm gcd */
{
  long unsigned int a=labs(a1);
  long unsigned int b=labs(b1);
  
  long unsigned int c;
  while(a!=0)
    {
      c=a;
      a=b%a;
      b=c;
    };
  return((int64_t) b);
};


// u,v,w <- sub determinants of th
void calc_uvw(int64_t *th, int64_t *uvw)
{
  
  printf("In calc_uvw with th = \n");
  for(uint64_t r=0;r<2;r++)
    {
      printf("[");
      for(uint64_t c=0;c<3;c++)
	printf(" %ld",th[3*r+c]);
      printf(" ]\n");
    }
  
  uvw[0]=th[1]*th[5]-th[2]*th[4];
  uvw[1]=th[0]*th[5]-th[2]*th[3];
  uvw[2]=th[0]*th[4]-th[1]*th[3];
  int64_t g1=gcd(labs(uvw[0]),labs(uvw[1]));
  int64_t g2=gcd(g1,labs(uvw[2]));
  uvw[0]/=g2;
  uvw[1]/=g2;
  uvw[2]/=g2;
  
}

void calc_YZ(int64_t *uvw, double *YZ)
{
  int64_t u,v,w;
  u=uvw[0];
  v=uvw[1];
  w=uvw[2];
  if(w>=0)
    {
      YZ[0]=w+v-u;      
      YZ[1]=(double)w/2.0+v-u;
    }
  else
    {
      YZ[0]=(double) w/2.0+v-u;
      YZ[1]=w+v-u;
    }
}

//A = [1,0,0]
//    [1,1,1]
//    [2,0,2]
void do_A(int64_t *th, int64_t *th1)
{
  th1[0]=th[0]+th[1]+2*th[2];
  th1[1]=th[1];
  th1[2]=th[1]+2*th[2];
  th1[3]=th[3]+th[4]+2*th[5];
  th1[4]=th[4];
  th1[5]=th[4]+2*th[5];
}

// BA = [0,1,0]
//      [2,0,1]
//      [2,0,2]
void do_BA(int64_t *th, int64_t *th1)
{
  th1[0]=th[1]*2+th[2]*2;
  th1[1]=th[0];
  th1[2]=th[1]+2*th[2];
  
  th1[3]=th[4]*2+th[5]*2;
  th1[4]=th[3];
  th1[5]=th[4]+2*th[5];
}

int main()
{
  int64_t *th,*th1;
  th=(int64_t *)malloc(sizeof(int64_t)*6);
  th1=(int64_t *)malloc(sizeof(int64_t)*6);
  th[0]=1;th[1]=1;th[2]=0;th[3]=0;th[4]=0;th[5]=1;
  
  for(;;)
    {
      int64_t uvw[3];
      calc_uvw(th,uvw);
      printf("u:%ld v:%ld w:%ld\n",uvw[0],uvw[1],uvw[2]);
      double YZ[2];
      calc_YZ(uvw,YZ);
      printf("Y:%f Z:%f\n",YZ[0],YZ[1]); 
      
      double Z=YZ[1];
      if(Z>=0.0) // do A
	{
	  do_A(th,th1);
	  printf("A\n");
	  int64_t *tmp=th;
	  th=th1;
	  th1=tmp;
	  continue;
	}

      double Y=YZ[0];
      if(Y<=0)
	{
	  do_BA(th,th1);
	  printf("BA\n");
	  int64_t *tmp=th;
	  th=th1;
	  th1=tmp;
	  continue;
	}
      printf("\nIndeterminate found.\n");
      break;
    }

  int64_t den=th[1]+2*th[2];
  int64_t num=th[4]+2*th[5];
  
  printf("\n%lu/%lu\n -> %14.12f\n",den,num,(double) den/ (double) num);
  return 0;
}
