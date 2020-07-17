/*
  File f_copy.cpp
  reads 1/2 FFT data from f_hat_even1.6
  and copies the complex conjugate halves
  suitable for input to f_hat_split1.2
*/

#include "../includes/int_double12.0.h"
#include "../includes/fft_defs.h"
#include "../includes/im_s.h"
#include "f_defs.h"
#include "stdio.h"

#define muassert(statement) if(!(statement)) {printf("Error in myassert. Exiting.\n");exit(0);}

void print_usage (const char *command )
{
  printf("Usage:- %s <file_format> <n files from f_hat_>. <1 even 0 odd>\n",command);
  exit(0);
}

int main(int argc, char **argv)
{
  _fpu_rndd();
  if(argc!=4)
    print_usage(argv[0]);
  int num_files=atoi(argv[2]);
  char buff[1024];
  FILE *infile,*outfile;
  int fn;
  unsigned int N,q,num_prims,num_ss,n0,index,num_ss1;
  int_complex omega,oz,*zs;
  long unsigned int lN;
  double eta;
  bool real_p,even_p,neg_one;
  even_p=(atoi(argv[3])!=0);
  for(fn=0;fn<num_files;fn++)
    {
      sprintf(buff,argv[1],fn); //
      printf("Processing file %s\n",buff);
      infile=fopen(buff,"rb");
      if(!infile)
	{
	  printf("Error opening file %s for binary input. Exiting.\n",buff);
	  exit(0);
	}
      sprintf(buff,argv[1],num_files*2-1-fn);
      outfile=fopen(buff,"wb");
      if(!infile)
	{
	  printf("Error opening file %s for binary input. Exiting.\n",buff);
	  exit(0);
	}
      muassert(fread(&q,sizeof(unsigned int),1,infile)==1);
      printf("q=%d\n",q);
      muassert(fwrite(&q,sizeof(unsigned int),1,outfile)==1);
      muassert(fread(&num_prims,sizeof(unsigned int),1,infile)==1);
      printf("num_prims=%d\n",num_prims);
      muassert(fwrite(&num_prims,sizeof(unsigned int),1,outfile)==1);
      muassert(fread(&N,sizeof(unsigned int),1,infile)==1);
      muassert(fwrite(&N,sizeof(unsigned int),1,outfile)==1);
      muassert(fread(&eta,sizeof(double),1,infile)==1);
      muassert(fwrite(&eta,sizeof(double),1,outfile)==1);
      muassert(fread(&num_ss,sizeof(unsigned int),1,infile)==1);
      if(fn==0)
	muassert(zs=(int_complex *)_aligned_malloc(sizeof(int_complex)*(num_ss+1),16));
      if((fn==0)||(fn==(num_files-1))) // 0'th file contains real f[0]
	                               // n-1 contains real f[n/2]
	num_ss1=num_ss-1;
      else
	num_ss1=num_ss;
      muassert(fwrite(&num_ss1,sizeof(unsigned int),1,outfile)==1);
      muassert(fread(&n0,sizeof(unsigned int),1,infile)==1);
      muassert(fwrite(&n0,sizeof(unsigned int),1,outfile)==1); // this is garbage
      for(unsigned int prim=0;prim<num_prims;prim++)
	{
	  //printf("prim=%d\n",prim);
	  // neg_one
	  muassert(fread(&neg_one,sizeof(bool),1,infile)==1);
	  muassert(fwrite(&neg_one,sizeof(bool),1,outfile)==1);
	  // real_p
	  muassert(fread(&real_p,sizeof(bool),1,infile)==1);
	  muassert(fwrite(&real_p,sizeof(bool),1,outfile)==1);
	  if(!real_p) prim++;
	  if(neg_one==even_p) continue;
	  // index
	  muassert(fread(&index,sizeof(unsigned int),1,infile)==1);
	  muassert(fwrite(&index,sizeof(unsigned int),1,outfile)==1);
	  // omega
	  muassert(fread(&omega,sizeof(int_complex),1,infile)==1);
	  muassert(fwrite(&omega,sizeof(int_complex),1,outfile)==1);
	  // zs
	  muassert(fread(zs,sizeof(int_complex),num_ss,infile)==num_ss);
	  if(fn==0)
	    {
	      for(unsigned int ptr=num_ss-1;ptr>0;ptr--) // don't write first element
		{
		  oz=conj(zs[ptr]);
		  muassert(fwrite(&oz,sizeof(int_complex),1,outfile)==1);
		}
	    }
	  else
	    {
	      if(fn==num_files-1) // dont write last element
		{
		  for(int ptr=num_ss-2;ptr>=0;ptr--)
		    {
		      oz=conj(zs[ptr]);
		      muassert(fwrite(&oz,sizeof(int_complex),1,outfile)==1);
		    }
		}
	      else
		for(int ptr=num_ss-1;ptr>=0;ptr--)
		  {
		    oz=conj(zs[ptr]);
		    muassert(fwrite(&oz,sizeof(int_complex),1,outfile)==1);
		  }
	    }
	  if(real_p) continue;
	  // index
	  muassert(fread(&index,sizeof(unsigned int),1,infile)==1);
	  muassert(fwrite(&index,sizeof(unsigned int),1,outfile)==1);
	  // omega
	  muassert(fread(&omega,sizeof(int_complex),1,infile)==1);
	  muassert(fwrite(&omega,sizeof(int_complex),1,outfile)==1);
	  // zs
	  muassert(fread(zs,sizeof(int_complex),num_ss,infile)==num_ss);
	  if(fn==0)
	    {
	      for(unsigned int ptr=num_ss-1;ptr>0;ptr--) // don't write first element
		{
		  oz=conj(zs[ptr]);
		  muassert(fwrite(&oz,sizeof(int_complex),1,outfile)==1);
		}
	    }
	  else
	    {
	      if(fn==num_files-1) // dont write last element
		{
		  for(int ptr=num_ss-2;ptr>=0;ptr--)
		    {
		      oz=conj(zs[ptr]);
		      muassert(fwrite(&oz,sizeof(int_complex),1,outfile)==1);
		    }
		}
	      else
		for(int ptr=num_ss-1;ptr>=0;ptr--)
		  {
		    oz=conj(zs[ptr]);
		    muassert(fwrite(&oz,sizeof(int_complex),1,outfile)==1);
		  }
	    }
	}
      fclose(infile);fclose(outfile);
    }
  return(0);
}
