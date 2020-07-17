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


uint64_t ptr=0,last_pix=0;
int_double low_sum;    
void callback(uint64_t p)
{
  while(ptr<p)
    low_pix[ptr++]=last_pix;
  last_pix++;
  low_sum*=(p-d_one)/p;
}

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
	  int_double x=int_double(x1,x2);
	  int_double pix=pistar(x);
	  pix+=int_double(ipix1,ipix2);
	  pix*=(x2-x1);
	  pix/=sqr(x);
	  res+=pix;
	  ptr++;
	  x1=x2;ipix1=ipix2;
	}
      printf("Integrated to x=%lu pi(x)=%lu ptr=%lu\n",x2,ipix2,ptr);
      print_int_double_str("Int=",res);
      sofar=low_sum-res+log(log(int_double(x2)))-(pistar(x2)+ipix2)/x2+d_one/(x2-1);
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

  low_sum=1.0;

  printf("Sieving small primes.\n");
  primesieve_callback_primes(2,SIEVE_LEN,callback);
  while(ptr<=SIEVE_LEN)
    low_pix[ptr++]=last_pix;

  low_sum=log(low_sum)+d_gamma-d_one/(X0-1)+(pistar(X0)+low_pix[X0])/X0;
  print_int_double_str("low sum = ",low_sum);


  printf("Integrating.\n");
  int_pistar(low_sum);

  return(0);

}

  
