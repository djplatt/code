/*
 *
 * Count how many triples of the form {a,b,d} with a<b<d, b^6<=d<=b^7.7
 * can be extended to {a,b,c,d} with b<c<d
 * We expect the answer to be zero
 */
#include "inttypes.h"
#include "stdio.h"
#include "stdlib.h"
#include "gmp.h"

mpz_t tmp1,tmp2,tmp3,ab1,r,t,s,d,upb,lwb;

mpz_t lwbc,upbc,tc,sc,rc,c;

// Given a,b and d, find all c in (b,d) such that {a,b,c,d} is a D. quadruple
uint64_t abd(uint64_t a,uint64_t b, int64_t ss)
{
  uint64_t nsols=0;
  uint64_t lwbc=b+1;
  mpz_mul_ui(tmp1,d,b);
  mpz_add_ui(tmp1,tmp1,1);
  mpz_sqrt(rc,tmp1); // rc=sqrt(bd+1)
  mpz_sub_ui(upbc,d,1);
  mpz_set_ui(tc,1); // starting column vector (tc,sc) is (1,1) or (1,-1)
  mpz_set_si(sc,ss);
  while(1)
    {
      mpz_mul(tmp1,tc,rc);
      mpz_mul(tmp2,sc,d);
      mpz_add(tmp3,tmp1,tmp2);
      mpz_mul_ui(tmp1,tc,b);
      mpz_mul(tmp2,rc,sc);
      mpz_add(sc,tmp1,tmp2); // sc=tc*b+sc*rc
      mpz_set(tc,tmp3); // tc=tc*rc+sc*d
      mpz_mul(tmp1,sc,sc);
      mpz_sub_ui(tmp1,tmp1,1);
      mpz_div_ui(c,tmp1,b); // c=(sc^2-1)/b=(tc^2-1)/d
      if(mpz_cmp_ui(c,lwbc)>=0) // c big enough
	{
	  if(mpz_cmp(c,upbc)<=0) // c not too big
	    {
	      mpz_mul_ui(tmp1,c,a);
	      mpz_add_ui(tmp1,tmp1,1);
	      if(mpz_perfect_square_p(tmp1)) // check ac+1 is a square
		{
		  //printf(" works\n");
		  nsols++;
		}
	      //else
	      //printf(" doesn't work\n");
	    }
	  else // c too big
	    break;
	}
    }
  return(nsols);
}

// given (a,b) such that ab+1=r^2, find all d such that {a,b,d} is a triple
// with d in [b^6,b^7.7]

//abp(a,b,r)={sols=0;ld=b^6;ud=floor(b^7.7);t=1;s=1;d=0;while(1,t1=t*r+s*b;s=a*t+r*s;t=t1;d=(s*s-1)/a;if(d>=ld,if(d<=ud,print(a," ",b," ",d);sols+=abd(a,b,d,floor(sqrt(b*d+1))),break)));return(sols)};
uint64_t abp(uint64_t a,uint64_t b,int64_t ss)
{
  uint64_t nsols=0;
  mpz_set_ui(t,1);
  mpz_set_si(s,ss); // vector [t,s]=[1,1] or [1,-1]
  mpz_ui_pow_ui(lwb,b,5);
  mpz_ui_pow_ui(upb,b,8);
  //mpz_ui_pow_ui(tmp2,b,77); // might be quicker to do ^8 instead
  //mpz_root(upb,tmp2,10);
  while(1)
    {
      mpz_mul(tmp1,t,r);
      mpz_mul_ui(tmp2,s,b);
      mpz_add(tmp3,tmp1,tmp2);
      mpz_mul_ui(tmp1,t,a);
      mpz_mul(tmp2,r,s);
      mpz_add(s,tmp1,tmp2);
      mpz_set(t,tmp3);
      mpz_mul(tmp1,s,s);
      mpz_sub_ui(tmp1,tmp1,1);
      mpz_div_ui(d,tmp1,a); //d=(s^2-1)/a=(t^2-1)/b
      if(mpz_cmp(d,lwb)>=0)
	{
	  if(mpz_cmp(d,upb)<=0) // we have {a,b,d} with d in range
	    {
	      //printf("%lu %lu ",a,b);
	      //mpz_out_str(NULL,10,d);
	      uint64_t dsols=abd(a,b,1)+abd(a,b,-1); // try to extend to {a,b,c,d}
	      //printf(" resulted in %lu solutions.\n",dsols);
	      nsols+=dsols;
	    }
	  else
	    break; // d too big, so break
	}
    }
  return(nsols);
}

uint64_t ab(uint64_t a, uint64_t b)
{
  mpz_sqrt(r,ab1);
  return(abp(a,b,1)+abp(a,b,-1));
}

      
int main(int argc, char **argv)
{
  if(argc!=4)
    {
      printf("Usage:- %s <a0> <a1> <b1>.\n",argv[0]);
      exit(0);
    }

  int64_t a0=atol(argv[1]);
  int64_t a1=atol(argv[2]);
  int64_t b1=atol(argv[3]);;

  //mpz_inits(tmp1,tmp2,tmp3,ab1,r,t,s,d,upb,lwb,NULL);
  
  mpz_init(tmp1);
  mpz_init(tmp2);
  mpz_init(tmp3);
  mpz_init(ab1);
  mpz_init(r);
  mpz_init(s);
  mpz_init(t);
  mpz_init(d);
  mpz_init(upb);
  mpz_init(lwb);
  mpz_init(upbc);
  mpz_init(lwbc);
  mpz_init(rc);
  mpz_init(sc);
  mpz_init(tc);
  mpz_init(c);  

  uint64_t a,b,nsols=0;
  for(a=a0;a<=a1;a++)
    for(b=a+3;b<=b1;b++)
      {
	mpz_set_ui(tmp1,a);
	mpz_mul_ui(ab1,tmp1,b);
	mpz_add_ui(ab1,ab1,1);
	if(mpz_perfect_square_p(ab1))
	  nsols+=ab(a,b);
      }
  printf("We found %lu solutions with a in [%ld,%ld], b<=%lu.\n",nsols,a0,a1,b1);
  return(0);
}
