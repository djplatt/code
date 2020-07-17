#include "stdio.h"
#include "stdlib.h"
#include "math.h"

long unsigned int Q;
double MAX_LOG,MIN_LOG;
long unsigned int NUM_BOXES;

long unsigned int *boxes;

int ignored=0;

double xdelta ()
{
  return((MAX_LOG-MIN_LOG)/NUM_BOXES);
}

double ydelta()
{
  return(2*M_PI/NUM_BOXES);
}

void put_box(double x, double y)
{
  long int xbox=(x-MIN_LOG)*NUM_BOXES/(MAX_LOG-MIN_LOG);
  long int ybox=(y+M_PI)*NUM_BOXES/(2*M_PI);
  //printf("%f going into box %ld\n",x,box);
  if((xbox<0)||(xbox>=NUM_BOXES)||(ybox<0)||(ybox>=NUM_BOXES))
    {
      printf("Ignoring %f %f \n",x,y);
      ignored++;
    }
  else
    boxes[xbox*NUM_BOXES+ybox]++;
}


inline double relog(double *z)
{
  double res=0.5*log((z[0]*z[0]+z[0]*z[0])/Q);
  //printf("relog(%f,%f)=%f\n",z[0],z[1],res);
  return(res);
}

inline double imlog(double *z)
{
  return(atan2(z[1],z[0]));
}

int main(int argc, char ** argv)
{
  double z[2],x,y,xdel,ydel;
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
  xdel=xdelta();
  ydel=ydelta();
  boxes=(long unsigned int *) malloc(sizeof(long unsigned int)*NUM_BOXES*NUM_BOXES);

  //fread(&Q,sizeof(long unsigned int),1,infile);
  Q=((long unsigned int) 1<<34)+25;
  for(i=0;i<NUM_BOXES*NUM_BOXES;i++)
    boxes[i]=0;

  fread(z,sizeof(double),2,infile); // principal so ignore it
  for(i=0;;i+=2)
    {
      if(!fread(z,sizeof(double),2,infile))
	break;
      else
	{
	  x=relog(z);
	  y=imlog(z);
	  put_box(relog(z),imlog(z));
	  if(z[1]!=0.0)
	    put_box(x,-y);
	}
    }
  // last one was real so
  printf("Last z was %f %f\n",z[0],z[1]);
  NUM_VALS=i-1;

  printf("I ignored %lu/%lu=%f%% of values.\n",ignored,NUM_VALS,(double)ignored/(double)NUM_VALS*100.0);
  for(i=0,x=MIN_LOG+xdel/2.0;i<NUM_BOXES;i++,x+=xdel)
    for(j=0,y=-M_PI+ydel/2.0;j<NUM_BOXES;j++,y+=ydel)
      printf("%f %f %10.8e\n",x,y,(double) boxes[i*NUM_BOXES+j]*NUM_BOXES*NUM_BOXES/NUM_VALS/(MAX_LOG-MIN_LOG)/2.0/M_PI);
  
  return(0);
}
