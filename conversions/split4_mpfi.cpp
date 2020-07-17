#include <stdio.h>
#include <stdlib.h>

#define in_num_s (512)
#define out_num_s (128)
#define fatal_error(st) {printf(st);printf("\n");exit(0);}
int main(int argc, char **argv)
{

  int i,j,num_s,file_N,rn,no_gaps,dbl_vec_len;
  double *dbl_vec;
  FILE *infile,*outfiles[4];
  if(argc!=6)
    fatal_error("split4_mpfi expects 5 arguments.");
  infile=fopen(argv[1],"rb");
  if(!infile)
    fatal_error("failed to open infile for binary input.");
  for(i=0;i<4;i++)
    {
      outfiles[i]=fopen(argv[2+i],"wb");
      if(!outfiles[i])
	fatal_error("failed to open outfile for binary output.");
    }


  fread(&num_s,sizeof(int),1,infile);
  if(num_s!=in_num_s)
    fatal_error("expected infile to be of length 512.");
  num_s>>=2;
  fread(&file_N,sizeof(int),1,infile);
  fread(&rn,sizeof(int),1,infile);
  if(rn!=1)
    fatal_error("expected rn to be 1.");
  fread(&no_gaps,sizeof(int),1,infile);
  dbl_vec_len=13+(no_gaps+1)*file_N*4;
  dbl_vec=(double *) malloc(sizeof(double)*dbl_vec_len);
  for(i=0;i<4;i++)
    {
      printf("Processing file number %d\n",i+1);
      fwrite(&num_s,sizeof(int),1,outfiles[i]);
      fwrite(&file_N,sizeof(int),1,outfiles[i]);
      fwrite(&rn,sizeof(int),1,outfiles[i]);
      fwrite(&no_gaps,sizeof(int),1,outfiles[i]);
      for(j=0;j<num_s;j++)
	{
	  fread(dbl_vec,sizeof(double),dbl_vec_len,infile);
	  fwrite(dbl_vec,sizeof(double),dbl_vec_len,outfiles[i]);
	}
      fclose(outfiles[i]);
    }
  fclose(infile);
  return(0);
}


