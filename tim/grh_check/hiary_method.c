#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <flint/acb_dirichlet.h>

void check_lam(dirichlet_group_t G, dirichlet_char_t chi, uint64_t q, int64_t prec)
{
  acb_t ls[2];
  acb_init(ls[0]);acb_init(ls[1]);
  acb_t s;
  acb_init(s);
  arb_set_ui(acb_realref(s),1);
  acb_dirichlet_l_jet(ls[0],s,G,chi,0,2,prec);
  printf("L(1) = ");acb_printd(ls[0],30);
  printf("\nL'(1) = ");acb_printd(ls[1],30);
  acb_div(ls[0],ls[1],ls[0],prec); // L'/L (1)
  printf("\nL'/L(1) = ");acb_printd(ls[0],30);printf("\n");
  acb_clear(s);
  acb_clear(ls[0]);
  acb_clear(ls[1]);
}
  
bool first_real_primitive_character(dirichlet_char_t chi, dirichlet_group_t G, uint64_t q)
{
  if((q%8)!=0) // q not divisible by 8
    dirichlet_char_log(chi,G,q-1);
  else
    dirichlet_char_log(chi,G,q/2-1);
  return dirichlet_char_is_primitive(G,chi)&&dirichlet_char_is_real(G,chi);
}

// only call when 8|q
bool second_real_primitive_character(dirichlet_char_t chi, dirichlet_group_t G, uint64_t q)
{
  uint64_t res=q>>3;
  if((q%4)==1)
    dirichlet_char_log(chi,G,q/2-1);
  else
    dirichlet_char_log(chi,G,3*q/4-1);
    
  return dirichlet_char_is_primitive(G,chi)&&dirichlet_char_is_real(G,chi);
} 

int main(int argc, char **argv)
{
  printf("Command:-");
  for(uint64_t n=0;n<argc;n++)
    printf(" %s",argv[n]);
  printf("\n");
  if(argc!=3)
    {
      printf("Usage:- %s <q> <prec>.\n",argv[0]);
      return 0;
    }
  uint64_t q=atol(argv[1]);
  int64_t prec=atol(argv[2]);

  if((q%16)==0)
    {
      printf("There are no primitive real characters modulo 16n.\n");
      return 0;
    }
  
  dirichlet_group_t G;

  dirichlet_char_t chi;

  dirichlet_group_init(G,q);
  printf("There are %lu/%lu primitive characters mod %lu.\n",dirichlet_group_num_primitive(G),dirichlet_group_size(G),q);
  dirichlet_char_init(chi,G);
  /*
  while(dirichlet_char_next(chi,G)>=0)
    {
      printf("Exp: %lu : ",dirichlet_char_exp(G,chi));
      if(dirichlet_char_is_primitive(G,chi))
	printf("primitive : ");
      else
	printf("imprimitive : ");
      if(dirichlet_char_is_real(G,chi))
	printf("real\n");
      else
	printf("complex\n");
    }
	
  */
  
  if(!first_real_primitive_character(chi,G,q))
    {
      printf("Can't find primitive real character for q=%lu. Exiting.\n",q);
      return 0;
    }

  check_lam(G,chi,q,prec);
  
  if((q%8)==0)
    {
      if(!second_real_primitive_character(chi,G,q))
	{
	  printf("Can't find second real character for q=%lu. Exiting.\n",q);
	  return 0;
	}
      
      check_lam(G,chi,q,prec);
    }
  
  return 0;
  
}
    
  
