char *djp_fnames;
FILE *filenames;
FILE *zerofile=NULL;
#include "inttypes.h"
uint64_t djp_num_its;
uint64_t djp_it,djp_zn;
double djp_st[2];
int64_t djp_zs[2];
mpfr_t z,dz;

#define OP_ACC (101)

void djp_in_bytes(mpfr_ptr t, FILE *infile)
{
  uint64_t a;
  uint32_t b;
  uint8_t c;
  int res;

  if(fread(&a,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (a). Exiting.\n");
      exit(0);
    }
  //printf("a=%lu\n",a);
  if(fread(&b,sizeof(uint32_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (b). Exiting.\n");
      exit(0);
    }
  //printf("b=%u\n",b);
  if(fread(&c,sizeof(uint8_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (c). Exiting.\n");
      exit(0);
    }
  mpfr_set_ui(t,c,GMP_RNDN);
  mpfr_mul_2ui(t,t,32,GMP_RNDN);
  mpfr_add_ui(t,t,b,GMP_RNDN);
  mpfr_mul_2ui(t,t,64,GMP_RNDN);
  mpfr_add_ui(t,t,a,GMP_RNDN);
  mpfr_div_2ui(t,t,OP_ACC,GMP_RNDN);
}

int djp_next_file()
{
  if(zerofile)
    fclose(zerofile);
  if(!filenames)
    {
      mpfr_init(z);
      mpfr_init(dz);
      filenames=fopen(djp_fnames,"r");
    }
  char fname[1024];
  if(fscanf(filenames,"%s\n",fname)!=1)
    return -1;
  printf("Going to open %s.\n",fname);fflush(stdout);
  zerofile=fopen(fname,"r");
  fread(&djp_num_its,sizeof(uint64_t),1,zerofile);
  fread(djp_st,sizeof(double),2,zerofile);
  fread(djp_zs,sizeof(int64_t),2,zerofile);
  djp_it=0;
  djp_zn=djp_zs[0];
  mpfr_set_d(z,djp_st[0],GMP_RNDN);
  return 1;
}

int djp_next_zero(arb_t zz, int64_t PREC)
{
  if(!zerofile)
    djp_next_file();
  if(djp_zn==djp_zs[1])
    {
      djp_it++;
      if(djp_it==djp_num_its)
	djp_next_file();
      else
	{
	  fread(djp_st,sizeof(double),2,zerofile);
	  fread(djp_zs,sizeof(int64_t),2,zerofile);
	  djp_zn=djp_zs[0];
	  mpfr_set_d(z,djp_st[0],GMP_RNDN);
	}
    }
  djp_in_bytes(dz,zerofile);
  djp_zn++;
  mpfr_add(z,z,dz,GMP_RNDN);
  arb_set_interval_mpfr(zz,z,z,PREC);
  arb_add_error_2exp_si(zz,-OP_ACC-1);
  return 1;
}
  
