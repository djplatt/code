#include "../includes/int_double12.0.h"
#define INTER_H (7.0/32.0)
#define twoB_num (64)
#define twoB_den (5)
#define one_over_two_B ((double) ((double) twoB_den/(double) twoB_num))

//#define one_over_two_B (5.0/64.0)
#define INTER_N (20)

int_double d_inter_err;

#define M ((double) 2.5)
#define ZETA_M (int_double(1.342,1.343))
#define ERR_CNST (int_double(2.378,2.379)) 

int_double G(int n, const double &t0, const int_double &B)
{
  int_double res=pow(int_double(1.5+t0+int_double((n+INTER_N)/2.0)/B),9.0/16.0);
  //print_int_double_str("res=",res);
  res*=exp(-sqr(int_double(n+INTER_N)/B/INTER_H)/8.0)/d_pi/(n+INTER_N);
  //printf("G(%d,%10.8e) returning",n,t0);
  //print_int_double_str("",res);
  return(res);
}

void set_d_inter_err(int q, double t0)
{
  int_double err_cnst=int_double(42.761,42.762);// =sqrt(Pi)*Zeta(9/8)*exp(1/6)*2^1.25
  int_double zeta_3=int_double(1.2020,1.2021);
  //printf("Setting d_inter_err for q=%d and t0=%10.8e\n",q,t0);
  int_double B=int_double(0.5)/one_over_two_B;
  int_double G0=G(0,t0,B);
  //print_int_double_str("G0=",G0);
  d_inter_err=pow(int_double(q)/d_pi,M/2.0)*2.0*zeta_3*exp(sqr(int_double(M)/INTER_H)/2-d_pi*2.0*M*B);
  d_inter_err*=INTER_H;
  d_inter_err*=(t0+int_double(INTER_H)/sqrt(d_two_pi)+1.0+int_double(0.5)/sqrt(int_double(2)))/M;
  //print_int_double_str("Ierr=",d_inter_err);
  int_double Eerr=G0/(d_one-G(1,t0,B)/G0)*pow(int_double(q)/d_two_pi,5.0/16.0)*err_cnst;
  //print_int_double_str("Eerr=",Eerr);
  d_inter_err+=Eerr;
  d_inter_err.left=d_inter_err.right;
  //print_int_double_str("inter_err set to",d_inter_err);
  //exit(0);
}

int main()
{
  _fpu_rndd();
  set_d_inter_err(3,1e8/3.0);
  print_int_double_str("",d_inter_err);
  set_d_inter_err(400000,1e8/400000.0);
  print_int_double_str("",d_inter_err);
  return(0);
}
