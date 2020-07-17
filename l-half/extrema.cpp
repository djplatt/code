#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define IGNORE_ME ((double) 1e10)
#define MAX_X ((double) 2.0)
#define MAX_Y ((double) 2.0)
#define NUM_DIVS (10)
#define X_DIV ((double) 2.0*MAX_X/NUM_DIVS)
#define Y_DIV ((double) 2.0*MAX_Y/NUM_DIVS)
#define MAX_VALS (1000000)
#define NUM_VALS (137438953680/16)

int boxes[(NUM_DIVS)*(NUM_DIVS)];

double in[MAX_VALS*2];

int ignored=0;

inline void put_box(double x, double y)
{
  int xpos,ypos;
  double x1,y1;
  printf("value=%f %f\n",x,y);
  return;
  if((x<=-MAX_X)||(x>=MAX_X)||(y<=-MAX_Y)||(y>=MAX_Y))
    {
      //printf("ignoring %f %f\n",x,y);
      ignored++;
      return;
    }
  x1=x+MAX_X;
  y1=y+MAX_Y;
  xpos=floor(x1/X_DIV);
  ypos=floor(y1/Y_DIV);
  //printf("Putting %f %f into box %ld %ld.\n",x,y,xpos,ypos);
  boxes[xpos*NUM_DIVS+ypos]++;
}

int main()
{
  double x,y;
  int i,j;
  float w,w1;
  FILE *r_infile=fopen("L_half_data.dat","rb");
  fread(in,sizeof(double),MAX_VALS,r_infile);
  fseek(r_infile,-(MAX_VALS+1)*sizeof(double),SEEK_END);
  fread(&in[MAX_VALS],sizeof(double),MAX_VALS,r_infile);
  fclose(r_infile);
  FILE *ofile=fopen("gnuplot.dat","wb");
  put_box(in[0],0.0);
  for(i=1;i<MAX_VALS;i++)
    {
      put_box(in[i],in[2*MAX_VALS-i-1]);
      put_box(in[i],-in[2*MAX_VALS-i-1]);
    }
  for(i=0;i<NUM_DIVS;i++)
    {
      for(j=0;j<NUM_DIVS;j++)
	printf("%7d ",boxes[i*NUM_DIVS+j]);
      printf("\n");
    }
  w=NUM_DIVS;
  fwrite(&w,sizeof(float),1,ofile);
  for(i=0,w=-MAX_Y+MAX_Y/NUM_DIVS;i<NUM_DIVS;i++,w+=MAX_Y*2.0/NUM_DIVS)
    fwrite(&w,sizeof(float),1,ofile);
  for(i=0,w=-MAX_X+MAX_X/NUM_DIVS;i<NUM_DIVS;i++,w+=MAX_X*2.0/NUM_DIVS)
    {
      fwrite(&w,sizeof(float),1,ofile);
      for(j=0;j<NUM_DIVS;j++)
	{
	  w1=boxes[i*NUM_DIVS+j];
	  fwrite(&w1,sizeof(float),1,ofile);
	}
    }
    
  printf("I ignored %f%% of values.\n",(double)ignored/(double)MAX_VALS*100.0);
  return(0);
}
