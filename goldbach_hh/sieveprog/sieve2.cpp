#include <stdlib.h>
#include <stdio.h>
#include <gmpxx.h>
#include <math.h>
#include "diophappr.h"
#include "solvemodz2.h"
#include "simplesieve.h"
/*  Compile using
g++ sieve.cpp -osieve -O2 -frounding-math -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lgmp */

void NewSegSiev(short *s, long n, long D, long K, long quot)
/*ensure: s[j] = (n-D+j is prime) for 0<=j<=2 D */
/* parameter kappa = 1/quot */
{
  mpz_class sqt, M2kap, sqtR;
  mpz_class npD(n+D);
  long M, R, m, m0, x; /*, trupos, falspos;*/
  long np;
  rat alpha0, alpha1;
  mpq_class eta;
  intlist *boldr, *r;
  
  mpz_sqrt(sqt.get_mpz_t(), npD.get_mpz_t());
  x = sqt.get_si();                     /* x = (int) sqrt(n+D) */
  		  
  SubSegSiev(s,n-D,2*D,K*D);
  /*  fprintf(stderr,"dudorino\n");*/
  for(M = K*D+1; M <=x; M+=(2*R+1)) {
     M2kap = ((((((mpz_class) M)*M)*D)/quot)/n);  /* (M^2)*D/(quot*n) */
    mpz_sqrt(sqtR.get_mpz_t(), M2kap.get_mpz_t());
    R = sqtR.get_si();  /* (int) (M*sqrt(kappa*D/n)) */
    
    m0 = M+R;
    alpha0.num = n%m0; alpha0.den = m0;
    alpha1.den = m0*m0; alpha1.num = (alpha1.den-n%alpha1.den);
    eta = (((mpq_class) D)/((mpq_class) M))*(1+1/((mpq_class) quot)); /* eta = 3D/2M */

    boldr = SolveModZ(alpha1,alpha0,R,eta);
    /*    trupos=falspos=0;*/
    /*    if(m0-R<=550481443 && m0+R>=550481443) {
      fprintf(stderr,"m0: %ld  R: %ld\n",M+R,R);
      fprintf(stderr,"%ld/%ld  %ld/%ld  %lg\n",
	      alpha1.num,alpha1.den,
	      alpha0.num,alpha0.den,
	      eta.get_d());
      fprintf(stderr,"%ld\n",(long) boldr);
      loud=1;
      boldr = SolveModZ(alpha1,alpha0,R,eta);
      loud=0;
      }*/
    for(r=boldr; r; r = r->next) {
      for(m = m0 + r->a; m<=m0+r->end; m+=r->q) {
      /*      if(m0-R<=550481443 && m0+R>=550481443)
	      fprintf(stderr,"%ld\t",r->it);
	      printf("m=%ld,\t",m); */
	np = ((n+D)/m)*m;
	if(np>=n-D && np<=n+D && np>m) /* { */ 
	  s[np-(n-D)]=0;
	  /*	printf("np = %ld, %lg\n",np,((double) (n%m))/m);
	  trupos++; 
	}  else {
	  	printf("nop = %lg\n",((double) (n%m))/m);
		falspos++;
	} */
      }
    }
    /*    fprintf(stderr,"true pos %ld, false pos %ld, neg %ld\n",trupos,falspos,2*R+1-trupos-falspos); */
  }
}

int main(int argc, char *argv[])
{
  long n, D, j;
  short *s;

  if(argc<=4) {
    if(argc<2) {
      n= 5000000000000000000;
      D=20000000;
    } else {
      n = atol(argv[1]);
      if(argc<3)
	D = 20000000;
      else
	D = atol(argv[2]);
    }

    s = (short *) calloc(2*D+1,sizeof(short));
    NewSegSiev(s,n,D,10,4);
  } else {
    /* call sieve2 400000000000 2000000 dummy
       (say) for control */
    n = atol(argv[1]);
    D = atol(argv[2]);
    SegSiev(s,n-D,2*D);
  }
  
  for(j=0; j<=2*D; j++)
    if(s[j])
      printf("%lu\n",n-D+j);
}
