#include "stdio.h"
#include "stdlib.h"
#include "math.h"

long unsigned int Q;
double MAX_LOG,MIN_LOG;
long unsigned int NUM_BOXES;

long unsigned int *boxes;

int ignored=0;

double delta ()
{
  return((MAX_LOG-MIN_LOG)/NUM_BOXES);
}

void put_box(double x)
{
  long int box=(long int) ((x-MIN_LOG)*(double) NUM_BOXES/(MAX_LOG-MIN_LOG));
  //printf("(x-MIN_LOG)*NUM_BOXES=%f\n",(x-MIN_LOG)*NUM_BOXES);
  //printf("%f going into box %ld\n",x,box);
  if((box<0)||(box>=NUM_BOXES))
    {
      //printf("Ignoring %f\n",x);
      ignored++;
    }
  else
    boxes[box]++;
}


inline double imlog(double *z)
{
  return(atan2(z[1],z[0]));
}

int main(int argc, char ** argv)
{
  double z[2],x,del;
  unsigned long int i,j,NUM_VALS;
  float w,w1;
  if(argc!=3)
    {
      printf("usage imlog <infile> <num boxes>\n");
      exit(0);
    }
  FILE *infile=fopen(argv[1],"rb");
  if(!infile)
    {
      printf("Error opening infile %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }
  MIN_LOG=-M_PI;
  MAX_LOG=M_PI;
  NUM_BOXES=atoi(argv[2]);
  del=delta();
  boxes=(long unsigned int *) malloc(sizeof(long unsigned int)*NUM_BOXES);

  fread(&Q,sizeof(long unsigned int),1,infile);
  //Q=((long unsigned int) 1<<34)+25;
  for(i=0;i<NUM_BOXES;i++)
    boxes[i]=0;

  fread(z,sizeof(double),2,infile); // principal so ignore it
  for(i=0;;i+=2)
    {
      if(!fread(z,sizeof(double),2,infile))
	break;
      put_box(imlog(z));
      if(z[1]==0.0)
	{
	  i--;
	}
      else
	{
	  z[1]=-z[1];
	  put_box(imlog(z));
	}
    }
  // last one was real so
  printf("Last z was %f %f\n",z[0],z[1]);
  NUM_VALS=i;

  printf("I ignored %lu/%lu=%f%% of values.\n",ignored,NUM_VALS,(double)ignored/(double)NUM_VALS*100.0);
  for(i=0,x=MIN_LOG+del/2.0;i<NUM_BOXES;i++,x+=del)
    printf("%f %10.8e\n",x,(double) boxes[i]*NUM_BOXES/NUM_VALS/(MAX_LOG-MIN_LOG));
  return(0);
}
