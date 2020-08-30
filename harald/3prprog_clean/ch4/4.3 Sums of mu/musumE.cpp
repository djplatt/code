/* Checks whether

|\sum_{n\leq x} \mu(n)/n|\leq sqrt(2/x)

using interval arithmetic.
*/

/* Compile with

g++ -fopenmp musumE.cpp -omusumE -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm

This assumes that crlibm is installed (in directory $CRDIR), that the
underlying processor is Intel, and that int_double14.2.h is in the
same directory as this file */

/* Test with

./musumE 1000000000 5 17 1

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "int_double14.2.h"
#include "mupsigma2.h"
#include "mubitarr.h"

double summu(ulong N, ulong m, int ntasks, ulong minsq, ulong minpr,
		 ulong outstep)
/* returns \sum_{n\leq N} \mu(n)/n */
/* stores the maximum of (\sum_{n\leq x} \mu(n)/n)*sqrt(x) for x\leq N+1 */
{
  ulong M, Mm, i, j,k,k2;
  ulong *imus, *isqf, *gcdind;
  unsigned char *pdiff, *indgcd;
  ulong **mus, **sqf;
  int th_id;
  double t0, t1;
  int_double psum, sum, *psumar;
  unsigned int i0;
  double maxrat, pmax, *pmaxar;
  long suma, *pasumar, pasuma;
  ulong *lor;
  
  sum = 0.0; suma = 0;

  initmubitdat(N,minsq,minpr,&imus, &isqf, &indgcd, &gcdind, &lor,
	       &pdiff, &M, &Mm);

  m = ((m-1)/M+1)*M;
  
  mus = (ulong **) calloc(ntasks,sizeof(ulong *));
  sqf = (ulong **) calloc(ntasks,sizeof(ulong *));
  psumar = (int_double *) calloc(ntasks,sizeof(int_double));
  pmaxar = (double *) calloc(ntasks,sizeof(double));
  pasumar = (long *) calloc(ntasks,sizeof(ulong));

  t0 = omp_get_wtime();

  maxrat=0.0;

  for(i=0; i<=N; i+=m*ntasks) {
 #pragma omp parallel private(k, k2, psum, pasuma, pmax, i0, th_id)
    {
#pragma omp for schedule(dynamic, 1)
      for(j=i; j<i+m*ntasks; j+=m) {
	pmax = 0.0; 
	i0 = ((j-i)/m);
	if(j<=N) {
	  mus[i0] = (ulong *) calloc(m/(8*sizeof(ulong))+1,sizeof(ulong));
	  sqf[i0] = (ulong *) calloc(m/(8*sizeof(ulong))+1,sizeof(ulong));

	  fillmubitblock(mus[i0], sqf[i0], imus, isqf, indgcd, gcdind, lor,
		      pdiff, j, m, minsq, minpr, M, Mm);
	  k2 = min(j+m,N+1);
	  psum = 0.0; pasuma = 0;
	  for(k=max(j,1); k<k2; k++) {
	    if(bitarrl(sqf[i0],k-j)) {
	      if(bitarrl(mus[i0],k-j)) {
		psum += 1.0/((int_double) k); pasuma++;
	      }
	      else {
		psum -= 1.0/((int_double) k); pasuma--;
	      }
	      if(max(-psum.left,-psum.right)>pmax)
		pmax = max(-psum.left,-psum.right);
	    }
	  }
	  psumar[i0] = psum; pmaxar[i0] = pmax; pasumar[i0] = pasuma;
       	  free(mus[i0]); free(sqf[i0]);
	} else {
	  psumar[i0] = 0.0; pmaxar[i0] = 0.0;
	}
      }
    }
 
    for(i0=0; i0<ntasks; i0++) {
      j = i + m*i0;
      if((((-((((int_double) max(-sum.right,-sum.left)) + pmaxar[i0]) *
	      sqrt((int_double) j+m)).right> maxrat)||(j+i0>N)) ||
	/* if there is any chance that an extremum is reached in the i0th
	   block, or if the i0th block is the last one... */	  
	  max(j,1)/outstep != (j+m-1)/outstep)
	 && j<=N)
	/* .. or if a multiple of outstep is contained in the i0th block
	   (we want to print out outstep and multiples of Mertens), then... */
	{
	  /*do exactly the same calculation as above */
	  
	  mus[0] = (ulong *) calloc(m/(8*sizeof(ulong))+1,sizeof(ulong));
	  sqf[0] = (ulong *) calloc(m/(8*sizeof(ulong))+1,sizeof(ulong));

	  fillmubitblock(mus[0], sqf[0], imus, isqf, indgcd, gcdind, lor,
		      pdiff, j, m, minsq, minpr, M, Mm);
			  
	  k2 = min(j+m,N+1);
	  for(k=max(j,1); k<k2; k++) {
	    if(bitarrl(sqf[0],k-j)) {
	      if(bitarrl(mus[0],k-j)) {
		sum += 1.0/((int_double) k); suma++;
	      } else {
		sum -= 1.0/((int_double) k); suma--;
	      }
	    }
	    if(-(max(-sum.left,-sum.right)*sqrt((int_double) k+1)).right
	       >maxrat && k>=3) {
	      maxrat = 
		-(max(-sum.left,-sum.right)*sqrt((int_double) k+1)).right;
	    }
	    if(!(k%outstep))
	      {
		printf("Mertens at %lu: %ld\n",k,suma);
		fflush(stdout);
	      }
	  }
	  free(mus[0]); free(sqf[0]);
      }
      else
	if(j<=N)
	  {sum += psumar[i0]; suma+=pasumar[i0];}
    }
  }
  
  t1 = omp_get_wtime();
  printf("Time spent: %g seconds\n",t1-t0);

  printf("m(%ld) lies in [%.15g %.15g]\n",N,sum.left,-sum.right);
  free(mus); freemubitdat(imus,isqf,indgcd,gcdind,lor,pdiff);
  free(psumar); free(pmaxar); free(pasumar); 

  return maxrat;
}

main(int argc, char *argv[])
{
  int nthreads, th_id;
  ulong N, m;
  ulong minsq, minpr;
  double maxr;

  // DJP Fix number of threads so we don't hoh LMFDB
  omp_set_num_threads(32);  
  //

   _fpu_rndd();
  
  if(argc>1)
    N = atol(argv[1]);
  else 
    N = 1000000000000;

  minsq=5; minpr = 17;
  
  if(argc>2)
    m = atol(argv[2]);
  else
    m = ((ulong) (-(4*sqrt((int_double) N)).right));

#pragma omp parallel private(th_id)
  {
    th_id = omp_get_thread_num();
    if ( th_id == 0 ) {
      nthreads = omp_get_num_threads();
      printf("There are %d threads\n",nthreads);
    }
  }

    
  maxr=summu(N, m, nthreads, minsq, minpr, 1000000000);
  printf("For 3<=x<=%ld, |m(x)| sqrt(x) is at most\n",N);
  printf("\t%.15g\n",maxr);
}


