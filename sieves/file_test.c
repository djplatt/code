#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

#define NUM_ITS ((long unsigned int) 2441407*8192)
#define BUFF_SIZE ((long unsigned int) 1)
#define fname "test.dat"
unsigned char buff[BUFF_SIZE];

int main()
{
  FILE *outfile,*infile;
  long unsigned int i;

  for(i=0;i<BUFF_SIZE;i++)
    buff[i]=i&255;

  printf("Going to write %lu bytes %lu times\n",BUFF_SIZE,NUM_ITS);

  time_t last_time=time(NULL);

  if(!(outfile=fopen(fname,"wb")))
    {
      printf("Failed to open file %s for output. Exiting.\n",fname);
      exit(0);
    }

  for(i=0;i<NUM_ITS;i++)
    fwrite(buff,sizeof(unsigned char),BUFF_SIZE,outfile);
  fclose(outfile);

  time_t this_time=time(NULL);
  printf("time to write file = %e\n",difftime(this_time,last_time));

  if(!(infile=fopen(fname,"rb")))
    {
      printf("Failed to open file %s for input. Exiting.\n",fname);
      exit(0);
    }
  for(i=0;i<NUM_ITS;i++)
    fread(buff,sizeof(unsigned char),BUFF_SIZE,infile);
  fclose(infile);

  last_time=time(NULL);
  printf("time to read file = %e\n",difftime(last_time,this_time));

  return(0);
}


  
