#include <stdio.h>
#include <stdlib.h>

#define fatal_error(st) {printf(st);printf("\n");exit(0);}

bool power_2(unsigned long int n)
{
  if(n==0)
    return(false);
  while(!(n&1))
    n>>=1;
  return(n==1);
}

int main(int argc, char **argv)
{

  int i,j,num_s,file_N,rn,no_gaps,dbl_vec_len,in_num_s;
  double *dbl_vec;
  FILE *infile,**outfiles;
  if(argc<2)
    fatal_error("Command line is <num files> <infile> <outfile1>..");
  int num_files=atoi(argv[1]);
  if(!power_2(num_files))
    fatal_error("Number of files must be a power of 2.");
  if(argc!=num_files+3)
    fatal_error("Command line is <num files> <infile> <outfile1>..");
  //printf("Splitting into %ld files.\n",num_files);
  infile=fopen(argv[2],"rb");
  if(!infile)
    fatal_error("failed to open infile for binary input.");
  outfiles=(FILE **)malloc(sizeof(FILE *)*num_files);
  for(i=0;i<num_files;i++)
    {
      outfiles[i]=fopen(argv[3+i],"wb");
      if(!outfiles[i])
	fatal_error("failed to open outfile for binary output.");
    }


  fread(&in_num_s,sizeof(int),1,infile);
  num_s=in_num_s/num_files;
  fread(&file_N,sizeof(int),1,infile);
  fread(&rn,sizeof(int),1,infile);
  if(rn!=1)
    fatal_error("expected rn to be 1.");
  fread(&no_gaps,sizeof(int),1,infile);
  dbl_vec_len=13+(no_gaps+1)*file_N*4;
  dbl_vec=(double *) malloc(sizeof(double)*dbl_vec_len);
  for(i=0;i<num_files;i++)
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


