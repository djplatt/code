#include "../includes/int_double12.0.h"
#include "inttypes.h"
#include "primesieve.h"

#define NFILES (9)
#define NRECS (900000)
#define SIEVE_LEN (1000000000L) // must be big enough to read sqrt(X)
#define X0 ((uint64_t) SIEVE_LEN) // this must live in filenames[0]
const char *filenames[NFILES]={
			 "/media/big_filesystem/data/pix/pix_kulsha_09.txt",
			 "/media/big_filesystem/data/pix/pix_kulsha_10.txt",
			 "/media/big_filesystem/data/pix/pix_kulsha_11.txt",
			 "/media/big_filesystem/data/pix/pix_kulsha_12.txt",
			 "/media/big_filesystem/data/pix/pix_kulsha_13.txt",
			 "/media/big_filesystem/data/pix/pix_kulsha_14.txt",
			 "/media/big_filesystem/data/pix/pix_kulsha_15.txt",
			 "/media/big_filesystem/data/pix/pix_kulsha_16.txt",
			 "/media/big_filesystem/data/pix/pix_kulsha_17.txt"};

uint64_t pixs[NFILES][NRECS+1];
uint64_t *low_pix;

/*
void read_pixs(uint64_t x0,uint64_t pix0)
{
  for(int i=0;i<NFILES;i++)
    {
      printf("Processing file %s.\n",filenames[i]);
      FILE* infile=fopen(filenames[i],"r");
      if(!infile)
	{
	  printf("Failed to open file %s. Exiting.\n",filenames[i]);
	  exit(0);
	}
      uint64_t x;
      if(fscanf(infile,"%lu %lu\n",&x,&pixs[i][0])!=2)
	{
	  printf("Error reading first line of file %s. Exiting.\n",filenames[i]);
	  exit(0);
	}
      if(x!=x0)
	{
	  printf("Expected file %s to start at x=%lu but read %lu. Exiting.\n",
		 filenames[i],x0,x);
	  exit(0);
	}
      if(pixs[i][0]!=pix0)
	{
	  printf("Expected file %s to start at pix=%lu but read %lu. Exiting.\n",
		 filenames[i],pix0,pixs[i][0]);
	  exit(0);
	}
      for(uint64_t j=1;j<=NRECS;j++)
	if(fscanf(infile,"%lu %lu\n",&x,&pixs[i][j])!=2)
	  {
	    printf("Error reading record %lu from file %s. Exiting.\n",
		   j,filenames[i]);
	    exit(0);
	  }
      x0=x;pix0=pixs[i][NRECS];
      fclose(infile);
    }
}
*/

uint64_t ptr=0,last_pix=0;
int_double low_sum;    
void callback(uint64_t p)
{
  while(ptr<p)
    low_pix[ptr++]=last_pix;
  last_pix++;
  low_sum+=log(1-d_one/p);
}
/*
int_double pix(uint64_t x)
{
  if(x<=PIX_START)
    return(int_double(low_pix[x],low_pix[x]));
  for(uint64_t i=0,y=PIX_START*10,z=10;i<NFILES;i++,y*=10,z*=10)
    {
      //printf("i=%lu y=%lu z=%lu\n",i,y,z);
      if(x<y)
	{
	  uint64_t ptr=(x-y/10)/z;
	  //printf("i=%lu y=%lu z=%lu ptr=%lu\n",i,y,z,ptr);
	  return(int_double(pixs[i][ptr],pixs[i][ptr+1]));
	}
    }
  printf("pix called with x=%lu. Too large! Exiting.\n",x);
  exit(0);
  return(0);
}
*/

int_double pistar(int_double x) // does everything but the pi(x)
{
  int_double res=0.0,lx=log(x);
  uint64_t n=2;
  while(true)
    {
      int_double mx=exp(lx/n);
      if(-mx.right<2.0)
	return(res);
      uint64_t l=floor(mx.left),r=floor(-mx.right);
      res+=int_double(low_pix[l],low_pix[r])/n;
      n++;
    }
}

// integrate from X0 to end of last file.
int_double int_pistar(int_double low_sum)
{
  int_double res=0.0,sofar;
  uint64_t file=0,ptr=0,x1,x2,ipix1,ipix2;
  while(file<NFILES)
    {
      FILE *infile=fopen(filenames[file],"r");
      x1=0;
      if(file==0)
	while(x1!=X0)
	  {
	    fscanf(infile,"%lu %lu\n",&x1,&ipix1);
	    ptr++;
	  }
      else
	{
	  fscanf(infile,"%lu %lu\n",&x1,&ipix1);
	  ptr++;
	}
      printf("Integrating from x=%li pi(x)=%lu ptr=%lu\n",x1,ipix1,ptr);
      while(ptr<=NRECS)
	{
	  fscanf(infile,"%lu %lu\n",&x2,&ipix2);
	  int_double pix=pistar(x);
	  pix+=int_double(ipix1,ipix2);
	  int_double x=int_double(x1,x2);
	  res+=(pix*(x2-x1))/sqr(x);
	  ptr++;
	  x1=x2;ipix1=ipix2;
	}
      printf("Integrated to x=%lu pi(x)=%lu ptr=%lu\n",x2,ipix2,ptr);
      print_int_double_str("Int=",res);
      sofar=low_sum-res+log(log(int_double(x2)))-(pistar(x2)+ipix2)/x2;
      print_int_double_str("So far=",sofar);
      fclose(infile);
      file++;
      ptr=0;
    }
  return(sofar);
}


int main(int argc, char**argv)
{
  printf("Command Line:- %s",argv[0]);
  for(int i=1;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");

  low_pix=(uint64_t *)malloc(sizeof(uint64_t)*(SIEVE_LEN+1));

  _fpu_rndd();

  low_sum=0.0;

  printf("Sieving small primes.\n");
  primesieve_callback_primes(2,SIEVE_LEN,callback);
  while(ptr<=SIEVE_LEN)
    low_pix[ptr++]=last_pix;

  low_sum+=d_gamma-d_one/(X0-1)+(pistar(X0)+low_pix[X0])/X0;
  print_int_double_str("low sum = ",low_sum);


  printf("Integrating.\n");
  int_pistar(low_sum);

  return(0);

}

  
