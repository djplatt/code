#include "stdio.h"
#include "stdlib.h"
#include "math.h"

long unsigned int Q;
double MAX_LOG,MIN_LOG,SQRT_Q;;
long unsigned int NUM_BOXES;

long unsigned int *boxes;

int ignored=0;

double delta ()
{
  return((MAX_LOG-MIN_LOG)/NUM_BOXES);
}

void put_box(double x)
{
  long int box=(x-MIN_LOG)*NUM_BOXES/(MAX_LOG-MIN_LOG);
  
  //printf("%f going into box %ld\n",x,box);
  if((box<0)||(box>=NUM_BOXES))
    {
      printf("Ignoring %f\n",x);
      ignored+=2;
    }
  else
    boxes[box]+=2;
}

void unput_box(double x)
{
  int box=(x-MIN_LOG)*NUM_BOXES/(MAX_LOG-MIN_LOG);
  if((box<0)||(box>=NUM_BOXES))
    ignored++;
  else
    boxes[box]--;
}

inline double absz(double *z)
{
  double res=sqrt(z[0]*z[0]+z[1]*z[1])*SQRT_Q;
  //printf("abs(%20.18e,%20.18e)=%f\n",z[0],z[1],res);
  return(res);
}

int main(int argc, char ** argv)
{
  double z[2],x,del;
  unsigned long int i,j,NUM_VALS;
  float w,w1;
  if(argc!=5)
    {
      printf("usage relog <infile> <min log> <max log> <num boxes>\n");
      exit(0);
    }
  FILE *infile=fopen(argv[1],"rb");
  if(!infile)
    {
      printf("Error opening infile %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }
  MIN_LOG=atof(argv[2]);
  MAX_LOG=atof(argv[3]);
  NUM_BOXES=atoi(argv[4]);
  del=delta();
  boxes=(long unsigned int *) malloc(sizeof(long unsigned int)*NUM_BOXES);

  fread(&Q,sizeof(long unsigned int),1,infile);
  //Q=((long unsigned int) 1<<34)+25;
  SQRT_Q=1.0/sqrt(Q);
  printf("Q=%lu\n",Q);
  for(i=0;i<NUM_BOXES;i++)
    boxes[i]=0;

  fread(z,sizeof(double),2,infile); // principal so ignore it
  for(i=0;;i+=2)
    {
      /*
      if((i%10000000)==0)
	printf("i=%lu\n",i);
      */
      if(!fread(z,sizeof(double),2,infile))
	break;
      else
	put_box(absz(z));
    }
  // last one was real so
  printf("Last z was %f %f\n",z[0],z[1]);
  unput_box(absz(z));
  NUM_VALS=i-1;

  printf("I ignored %lu/%lu=%f%% of values.\n",ignored,NUM_VALS,(double)ignored/(double)NUM_VALS*100.0);
  for(i=0,x=MIN_LOG+del/2.0;i<NUM_BOXES;i++,x+=del)
    printf("%f %10.8e\n",x,(double) boxes[i]*NUM_BOXES/NUM_VALS/(MAX_LOG-MIN_LOG));
  return(0);
}
