#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>
#include <gmp.h>
#include <mpfr.h>
#include <string.h>


#include "asm/config.h"
#include "if.h"

static char *zero_files_filename = NULL;
static FILE *zero_files_file = NULL;
static char *last_zfn=NULL;
static char *previous_last_zfn=NULL;


static char*
next_zero_filename(void)
{
  char zfname[201];

  for(;;) {
    if(feof(zero_files_file)) return NULL;
    if(fscanf(zero_files_file,"%200s\n",zfname)!=1)
      complain("Garbled line in %s\n",zero_files_filename);
    if(previous_last_zfn!=NULL) {
      if(strcmp(previous_last_zfn,zfname)==0) {
	free(previous_last_zfn);
	previous_last_zfn=NULL;
      }
      continue;
    }
    return strdup(zfname);
  }
}

/*
inline uint64_t reverse_me(uint64_t n)
{
  return __builtin_bswap64(n);
}

inline static u64_t
ru64big(FILE *f)
{
  u64_t r;
  if(fread(&r,sizeof(u64_t),1,f)!=1)
    complain("ru64big: premature eof, or I/O error\n");
  return reverse_me(r);
}
*/

static u64_t
ru64big(FILE*f)
{
  u64_t r,k;
  int c;
  for(r=0,k=0;k<8;k++) {
    c=getc(f);
    if(c<0) complain("ru64big: premature eof, or I/O error\n");
    r=(r<<8)+c;
  }
  return r;
}


static char *infn;
static FILE *inf=NULL;
static mpz_t range_start,range_end,checksum,zero128,d,mask;
static u64_t inf_nz=0;


static int
new_zeros_file(void)
{
  FILE *old_inf;

  if((old_inf=inf)!=NULL) {
    mpz_and(checksum,checksum,mask); // fix by djp since checksum can
                                       // overflow 128 bits

    if(mpz_cmp_ui(checksum,0)!=0)
      {
	fprintf(stderr,"WARNING: checksum error in %s\nchecksum was 0x",infn);
	mpz_out_str(stderr,16,checksum);printf("\n");
      }

    fclose(inf);
    inf=NULL;
    if(last_zfn!=NULL && strcmp(infn,last_zfn)==0) {
      free(infn);
      return -1;
    }
    free(infn);
  }

  if((infn=next_zero_filename())==NULL)
    return -1;
  if((inf=fopen(infn,"r"))==NULL) {
    perror(infn);
    complain("Could not open input file\n");
  }

  logbook(0,"%s ",infn);
  mpz_set_ui(range_start,ru64big(inf));
  mpz_mul_2exp(range_start,range_start,64);
  mpz_add_ui(range_start,range_start,ru64big(inf));
  if(old_inf!=NULL && mpz_cmp(range_start,range_end))
    fprintf(stderr,
	    "%s: WARNING: range_start differs from previous range_end\n",infn);

  mpz_set_ui(range_end,ru64big(inf));
  mpz_mul_2exp(range_end,range_end,64);
  mpz_add_ui(range_end,range_end,ru64big(inf));
  inf_nz=ru64big(inf);
  //printf("File %s contains %lu zeros.\n",infn,inf_nz);
  ru64big(inf);
  mpz_set_ui(checksum,ru64big(inf));
  mpz_mul_2exp(checksum,checksum,64);
  mpz_add_ui(checksum,checksum,ru64big(inf));
  mpz_set(zero128,range_start);
  return 1;
}


static int
read_zero(void)
{
  int c;

  while(inf_nz--==0)
      if(new_zeros_file()<0) return -1;

  if((c=getc(inf))<0)
    complain("%s: premature eof or I/O error\n",infn);
  if(c==0xFF) {
    if((c=getc(inf))<0)
      complain("%s: premature eof or I/O error\n",infn);
    mpz_set_ui(d,c);
    mpz_mul_2exp(d,d,64);
  } else {
    ungetc(c,inf);
    mpz_set_ui(d,0);
  }
  mpz_add_ui(d,d,ru64big(inf));
  mpz_add(zero128,zero128,d);
  mpz_sub(checksum,checksum,zero128);
  return 1;
}



int zs_get_next_zero(mpfr_t zero)
{
  if(read_zero()<0) return -1;
  mpfr_set_z(zero,zero128,GMP_RNDN);
  mpfr_div_2ui(zero,zero,64,GMP_RNDN);
  return 1;
}

int zs_init(char *filename)
{
	if(filename!=NULL){
		if((zero_files_file = fopen(filename, "r"))==NULL) return -1;
	}	
	else if((zero_files_file = fopen(zero_files_filename, "r"))==NULL) return -1;
	mpz_init(range_start);
	mpz_init(range_end);
	mpz_init(zero128);
	mpz_init(checksum);
	mpz_init(d);
	mpz_set_ui(range_end,0);
	mpz_init(mask);
	mpz_set_ui(mask,1);
	mpz_mul_2exp(mask,mask,128);
	mpz_sub_ui(mask,mask,1);
	//printf("mask=0x");mpz_out_str(stdout,16,mask);printf("\n");
	return 1;
}

int zs_close(void)
{
	if(zero_files_file!=NULL) fclose(zero_files_file);
	if(inf!=NULL) fclose(inf);
	inf=NULL;
	zero_files_file = NULL;
	free(last_zfn);
	last_zfn=NULL;
	free(previous_last_zfn);
	previous_last_zfn=NULL;
	inf_nz = 0;
	mpz_clear(range_start);
	mpz_clear(range_end);
	mpz_clear(zero128);
	mpz_clear(checksum);
	mpz_clear(d);
	return 0;
}

