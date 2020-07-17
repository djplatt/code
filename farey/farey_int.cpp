// Keith Briggs 2015-06-27

/*
def farey(n,asc=True):
  """Python function to print the nth Farey sequence, either ascending or descending."""
  if asc:
    a, b, c, d = 0, 1,  1  , n
  else:
    a, b, c, d = 1, 1, n-1 , n
  print("%d/%d" % (a,b),end=' ')
  while (asc and c <= n) or (not asc and a > 0):
    k = int((n + b)/d)
    a, b, c, d = c, d, k*c - a, k*d - b
    print("%d/%d" % (a,b),end=' ')
  print()
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <gmp.h>
#include "A005728_table.h" // from https://oeis.org/A005728/b005728.txt
#include "../includes/int_double13.0.h"

#define CPU_TIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + \
  ruse.ru_stime.tv_sec + 1e-6 * \
  (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))

typedef unsigned long ulong;
struct rusage ruse;
double t0,t1;

ulong A005728(ulong n) {
  // A005728 Number of fractions in Farey series of order n.
  // a(5)=11 because the fractions are 0/1, 1/5, 1/4, 1/3, 2/5, 1/2, 3/5, 2/3, 3/4, 4/5, 1/1.
  // a(n) = n(n+3)/2 - Sum(k = 2 to n, a([n/k])) = 1 + Sum_{i=1..n} phi(i).
  // a[0]=1, a[1]=2, 3, 5, 7, 11, 13, 19, 23, 29, 33, 43, 47, 59, 65, 73, 81,
  ulong k,sum;
  // replace below by global table from A005728_table.h
  //unsigned int A005728_table[]={1,2,3,5,7,11,13,19,23,29,33,43,47,59,65,73,81,97,103,121,129,141,151,173,181,201,213,231,243,271,279,309,325,345,361,385,397,433,451,475,491,531,543,585,605,629,651,697,713,755,775,807,831,883,901,941,965};
  if (n<A005728_table_size) return A005728_table[n];
  sum=0;
  for (k=2; k<=n; k++) sum+=A005728(n/k);
  return (n*(n+3))/2-sum;
}

void print_farey_discrepancy_ulong(ulong n) {
  // a,b,c,d <- c,d,k*c-a,k*d-b
  /*
     [a]'=[ 0  0       1       0] [a]
     [b] =[ 0  0       0       1] [b]
     [c] =[-1  0 (n+b)/d       0] [c]
     [d] =[ 0 -1       0 (n+b)/d] [d]
  */
  long l=0,a=0,b=1,c=1,d=n,cc,dd,k,m=A005728(n)-1;
  printf("%.2f seconds computing |F_n|\n",CPU_TIME-t0);
  int_double dist,dm=m,disc1=0.0,disc2=0.0;
  if (1==n) printf("# n\tdisc1\tdisc2\tdisc1/sqrt(n)\tdisc2*n\tcumulative time (seconds)\n");
  while (1) {
    cc=c; dd=d;
    k=(n+b)/d;
    d=k*d-b;
    c=k*c-a;
    a=cc; b=dd;
    //fprintf(stderr,"%lu/%lu ",a,b);
    l++;
    // the next line exploits symmetry around 1/2 to halve the runtime...
    if (2*l>=m) { disc1*=2.0; disc2*=2.0; break; } // 2*l==m fails for n=1
    //dist=a;dist/=b;dist-=(double)l/dm;
    cc=a*m-b*l;dist=cc;dist/=b*m;
    //printf("%10.8e %10.8e\n",dist,dist1);
    disc2+=dist*dist;      // Franel RH <=> O(n^{-1+epsilon})
    disc1+=abs(dist);     // Landau RH <=> O(n^{1/2+epsilon})
  }
  // columns 4 and 5 tend to constant <=> RH...
  t1=CPU_TIME;
  printf("%lu\t",n);print_int_double(disc1);printf("\t");
  print_int_double(disc2);printf("\t");
  print_int_double(disc1/sqrt(int_double(n)));printf("\t");
  print_int_double(disc2*n);
  printf("\t%.2f\n",t1-t0);
  fflush(stdout);
}
int main(int argc, char* argv[]) {
  ulong n,n0=1,n1=1000;
  if (argc>1) n0=atoi(argv[1]);
  if (argc>2) n1=atoi(argv[2]);
  _fpu_rndd();
  t0=CPU_TIME;
  for (n=n0; n<=n1; n++) {
    print_farey_discrepancy_ulong(n);
  }
  return 0;
}
