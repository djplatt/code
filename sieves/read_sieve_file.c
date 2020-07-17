#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "inttypes.h"
#include "../includes/pi_x.h"

#define debug printf("Reached line number %d.\n",__LINE__)


#define num_bits (sizeof(ptype)<<3)
#define log_num_bits (6)
#define all_ones ((num_bits<<1)-1)

#define primep(p,v) (v[p>>(log_num_bits+1)]&mask1[p&all_ones])
//#define primep(p,v) (v[p/(num_bits*2)]&mask1[p&all_ones])
#define clear_prime(p,v) (v[p>>(log_num_bits+1)]&=mask[p&all_ones])
//#define clear_prime(p,v) (v[p/(num_bits*2)]&=mask[p&all_ones])

inline void fatal_error(const char *str)
{
  fputs(str,stderr);
  fputs(" Exiting.\n",stderr);
  abort();
}

void print_usage()
{
  exit(0);

}


int main(int argc, char **argv)
{
  FILE *infile=fopen("sieve_546800.dat.old","rb");
  if(!infile)
    exit(0);
  FILE *outfile=fopen("sieve_547400.dat","wb");
  // will only contain data for iterations 547400-547492 inclusive
  ptype i,j,num_its,num_segs,count,del_sum,it,izero[2]={0,0};
  bigint start,del_sum2,bzero=0;
  fread(&num_its,sizeof(ptype),1,infile);
  fwrite(&num_its,sizeof(ptype),1,outfile);
  printf("Number of iterations=%lu\n",num_its);
  fread(&num_segs,sizeof(ptype),1,infile);
  fwrite(&num_segs,sizeof(ptype),1,outfile);
  printf("Number of segments per sieve=%lu\n",num_segs);
  
  for(i=0;i<num_its;i++)
    {
	fread(&it,sizeof(ptype),1,infile);
	fwrite(&it,sizeof(ptype),1,outfile);
	it=it-i+i*num_segs/2;
	fread(&start,sizeof(bigint),1,infile);
	fwrite(&start,sizeof(bigint),1,outfile);

	for(j=0;j<num_segs/2;j++)
	  {
	    //printf("starting at ");print_bigint(start);printf("\n");
	    start=start+((long unsigned int) 1<<32);
	    //printf("Iteration=%lu\n",it);
	    fread(&count,sizeof(ptype),1,infile);
	    //printf("count=%lu\n",count);
	    fread(&del_sum,sizeof(ptype),1,infile);
	    //printf("del_sum=%ld\n",del_sum);
	    fread(&del_sum2,sizeof(bigint),1,infile);
	    //printf("del_sum2=");print_bigint(del_sum2);printf("\n");
	    if((it>=547400)&&(it<=547492))
	      {
		printf("Writing iterations %ld starting at ",it);
		print_bigint(start);
		printf("\n");
		fwrite(&count,sizeof(ptype),1,outfile);
		fwrite(&del_sum,sizeof(ptype),1,outfile);
		fwrite(&del_sum2,sizeof(bigint),1,outfile);
	      }
	    else
	      {
		fwrite(izero,sizeof(ptype),2,outfile);
		fwrite(&bzero,sizeof(bigint),1,outfile);
	      }
	    it++;
	  }
	printf("\n");
    }



  return(0);
}
