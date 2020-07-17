#include "../includes/upsamdefs.h"

typedef complex<double> mcomplex;



bool testing=false;


void fatal_error (const char *str)
{
	printf("%s",str);
	exit(1);
}

mcomplex *in_array;
mcomplex *out_array;
fftw_plan p_fwd,p_bkwd;

// this contains the gaussian factors to use when upsampling
double gaussians[UP_SAMPLE_IN_WIDTH/2];


int_double expsincs[(UP_SAMPLE_RATE-1)*INTER_N*2];
int_double twoBpi;
// consists of UP_SAMPLE_RATE-1 "rows" of gauss*sinc
// first row centred at t0-n/2B = 1/UP_SAMPLE_RATE
// last row at                  = (UP_SAMPLE_RATE-1)/UP_SAMPLE_RATE
void setexpsincs()
{
	int offset,i,ptr=0;
	double del;
	int_double ex,si;

	for(offset=1;offset<UP_SAMPLE_RATE;offset++)
	{
		del=(((double) offset)/((double) UP_SAMPLE_RATE)+(double) INTER_N - 1.0)*one_over_two_B;
		for(i=0;i<INTER_N*2;i++,del-=one_over_two_B)
		{

			ex=exp(-int_double(del*del)/TWO_H_SQR);

			si=sincpi(twoBpi*del);

			expsincs[ptr++]=ex*si;
		}
	}
}

// please note double exp(double) does not work if you mess with rounding
// ergo use double exp1(double)
void set_gaussians()
{
	double i_over_h;
	for(int i=0;i<UP_SAMPLE_IN_WIDTH/2;i++)
	{
		i_over_h=(double) i/UP_SAMPLE_H;
		gaussians[i]=exp1(-(i_over_h*i_over_h));
	}
}

double exppittwos[UP_SAMPLE_SACRIFICE];
int_double exppittwos_d[UP_SAMPLE_SACRIFICE];

void setexppittwos(im_s *im_s_vec)
{
	for(int i=0;i<UP_SAMPLE_SACRIFICE;i++)
	{
	  //exppittwos[i]=exp1(-im_s_vec[i+1].im_s*d_pi.left/2);
	  //exppittwos_d[i]=exp(d_pi*im_s_vec[i+1].im_s/(-2.0));
	  exppittwos[i]=exp1(-(i+1)*one_over_two_B*d_pi.left/2);
	  exppittwos_d[i]=exp(d_pi*(i+1)*one_over_two_B/(-2.0));
	}
}

void setup(im_s* im_s_vec)
{
  twoBpi=d_one/one_over_two_B; // use sincpi from now on
  setexpsincs();
  setexppittwos(im_s_vec);
  __SSE_setcw(old__SSE_cw);
  set_gaussians();
  in_array=(mcomplex*) fftw_malloc(sizeof(fftw_complex) * UP_SAMPLE_IN_WIDTH);
  out_array=(mcomplex*) fftw_malloc(sizeof(fftw_complex) * UP_SAMPLE_OUT_WIDTH);
  
  p_fwd=fftw_plan_dft_1d(UP_SAMPLE_IN_WIDTH, reinterpret_cast<fftw_complex*>(in_array),
			 reinterpret_cast<fftw_complex*>(in_array), FFTW_FORWARD, FFTW_PLAN_MODE);
  
  p_bkwd=fftw_plan_dft_1d(UP_SAMPLE_OUT_WIDTH, reinterpret_cast<fftw_complex*>(out_array),
			  reinterpret_cast<fftw_complex*>(out_array), FFTW_BACKWARD, FFTW_PLAN_MODE);
  __SSE_setcw(new__SSE_cw);
}

void rig_setup(im_s* im_s_vec)
{
  twoBpi=d_one/one_over_two_B; // use sincpi
  setexpsincs();
  setexppittwos(im_s_vec);
}

// in_array[0..UP_SAMPLE_IN_WIDTH-2] contains the sampled data
// this is multiplied by the Gaussian
// the in_array[UP_SAMPLE_WIDTH-1] is set to avg of [0], [U_S_W-2]
// then FFT'd in place
// then zero padded into out_array
// then inverse FFT'd inplace
void up_sample ()
{

	int i,j;
	__SSE_setcw(old__SSE_cw);
	for(i=0;i<UP_SAMPLE_IN_WIDTH/2-1;i++)
	{
		in_array[UP_SAMPLE_IN_WIDTH/2-i-1]*=gaussians[i];
		in_array[UP_SAMPLE_IN_WIDTH/2+i]*=gaussians[i];
	}
	in_array[0]*=gaussians[i];
	in_array[UP_SAMPLE_IN_WIDTH-1]=(in_array[0]+in_array[UP_SAMPLE_IN_WIDTH-2])/2.0;

	fftw_execute(p_fwd);

	// put in first half of frequencies
	for(i=0;i<(UP_SAMPLE_IN_WIDTH>>1);i++)
		out_array[i]=in_array[i];
	// remember where we are in in_array
	j=i;
	// zero out central part of frequencies
	for(;i<(UP_SAMPLE_OUT_WIDTH-(UP_SAMPLE_IN_WIDTH>>1));i++)
		out_array[i]=mcomplex(0.0,0.0);
	// put in second half of frequencies
	for(;i<UP_SAMPLE_OUT_WIDTH;i++)
		out_array[i]=in_array[j++];

	fftw_execute(p_bkwd);
	__SSE_setcw(new__SSE_cw);
}

inline sign_t sign(int_double &x)
{
	if(x.left>0.0)
		return(POS);
	if(x.right>0.0)
		return(NEG);
	return(UNK);
}

inline sign_t d_sign(double t)
{
	if(t<0.0)
		return(NEG);
	if(t>0.0)
		return(POS);
	return(UNK);
}

void print_sign(sign_t sign)
{
	if(sign==POS)
		printf("+");
	else
	{
		if(sign==NEG)
			printf("-");
		else
			printf("?");
	}
}

void d_print_signs( int n, mcomplex *re_zs,  int step)
{
	int i;
	sign_t last_sign,this_sign;
	last_sign=d_sign(real(re_zs[0]));
	print_sign(last_sign);
	for(i=step;i<n;i+=step)
	{
		this_sign=(d_sign(real(re_zs[i])));
		if(this_sign!=last_sign)
		{
			//printf("\n");
			last_sign=this_sign;
		}
		print_sign(this_sign);
	}
	printf("\n");
}

void print_signs( int n, int_double *re_zs,  int step)
{
	int i;
	print_sign(sign(re_zs[0]));
	for(i=step;i<n;i+=step)
		print_sign(sign(re_zs[i]));
	printf("\n");
}

void print_zeros(int n, int_double *re_zs,im_s *im_s_vec)
{
	int i=1;
	sign_t this_sign,last_sign=sign(re_zs[0]);

	while(i<n)
	{
		this_sign=sign(re_zs[i]);
		if(this_sign==UNK)
		{i++;continue;}
		if(this_sign!=last_sign)
		{
		  printf("%7.4f\n",i*one_over_two_B);//im_s_vec[i].im_s);
			last_sign=this_sign;
		}
		i++;
	}
}

int check_count=0;

inline int_double check_sign1(int n0, int offset, int_double *re_zs1, int_double *re_zs2)
{
	int n,n1;
	int_double *expsinc,res=d_inter_err;
	check_count++;
	while(offset>UP_SAMPLE_RATE)
	{
		n0++;
		offset-=UP_SAMPLE_RATE;
	}
	while(offset<0)
	{
		n0--;
		offset+=UP_SAMPLE_RATE;
	}
	// t is now between [n0] and [n0+1]
	expsinc=&expsincs[(offset-1)*INTER_N*2];
	if(n0-INTER_N+1<0)
	{
		//if(sign(re_zs1[0])==sign(re_zs2[0]))
	  for(n=n0-INTER_N+1,n1=0;n<0;n++,n1++)
	    res+=re_zs2[-n]*expsinc[n1]*exppittwos_d[-n-1];
		//else
		//for(n=n0-INTER_N+1,n1=0;n<0;n++,n1++)
		//res+=(-re_zs2[-n])*expsinc[n1]*exppittwos_d[-n-1];
	  for(n=0;n<=n0+INTER_N;n++,n1++)
	    res+=re_zs1[n]*expsinc[n1];
		//print_int_double_str("check_sign res=",res);
		//sign_res=sign(res);
		//if(sign_res==UNK)
		//{
		//printf("In check_sign with n0=%d offset=%d\n",n0,offset);
		//print_int_double_str("Check sign calculation of unknown sign.",res);
		//}
	  return(res);
	}
	for(n=n0-INTER_N+1,n1=0;n<=n0+INTER_N;n++,n1++)
	  {
	    res+=re_zs1[n]*expsinc[n1];
	  }
	return(res);
}

inline sign_t check_sign(int n0, int offset, int_double *re_zs1, int_double *re_zs2)
{
	int_double res=check_sign1(n0,offset,re_zs1,re_zs2);
	//printf("In check_sign with n0=%d offset=%d\n",n0,offset);
	//print_int_double_str("Returning ",res);
	return(sign(res));
}
/*
// check for changes of sign between st and en inclusive
// detects ++++????+++++ and -----????---- as probable missed pairs
// outputs location of such regions for processing by upsamdouble
int rig_num_zeros1(int st, int en, int_double *re_zs1, int_double *re_zs2, FILE *outfile, q_state &qs)
{
  int count=0,ptr=st,offset=1,unk_ptr,unk_offset;
  sign_t last_sign,this_sign;
  bool unk_p=false;

  last_sign=sign(re_zs1[ptr]);
  if((last_sign==UNK)&&(st==0))
    {
      printf("F(0,chi) is of unknown sign.\n");
      print_int_double_str("F(0,chi)=",re_zs1[0]);
      print_int_double_str("F(1/2B,chi)=",re_zs1[1]);
      print_int_double_str("F(0,chi_bar)=",re_zs2[0]);
      print_int_double_str("F(1/2B,chi_bar)=",re_zs2[1]);
      //qs.type=RE_Z_0_PROB;
      //save_it(qs,outfile);
    }
  for(;last_sign==UNK;last_sign=check_sign(ptr,offset++,re_zs1,re_zs2));
  while(ptr<en)
    {
      if(offset==UP_SAMPLE_RATE)
	{
	  offset=1;
	  ptr++;
	  this_sign=sign(re_zs1[ptr]);
	}
      else
	this_sign=check_sign(ptr,offset++,re_zs1,re_zs2);
      if(this_sign==UNK)
	{
	  if(!unk_p)
	    {
	      unk_p=true;
	      unk_ptr=ptr;
	      unk_offset=offset-1;
	    }
	  continue;
	}
      if(this_sign==last_sign)
	{
	  if(unk_p)
	    {
	      qs.type=CHECK_PROB;
	      qs.n0=unk_ptr;qs.n02=ptr;
	      qs.offset1=unk_offset;qs.offset2=offset-1;
	      printf("Sign strayed into unknown and back near n0=%d-%d offset=%d-%d.",unk_ptr,ptr,unk_offset,offset-1);
	      if(this_sign==POS)
		{
		  printf(" Look for a negative.\n");
		  qs.exp_sign=NEG;
		}
	      else
		{
		  printf(" Look for a positive.\n");
		  qs.exp_sign=POS;
		}
	      save_it(qs,outfile);
	      unk_p=false;
	    }
	  continue;
	}
      unk_p=false;
      count++;
      last_sign=this_sign;
    }
  return(count);
}
*/

inline bool new_zero(sign_t *signs)
{
  int num_changes=0;
  sign_t last_sign=signs[0];
  //printf("In new_zero with ");
  //for(int i=0;i<5;i++)
  //print_sign(signs[i]);
  //printf("\n");
  for(int ptr=1;ptr<5;ptr++)
    {
      if(signs[ptr]==UNK)
	continue;
      if(signs[ptr]==last_sign)
	continue;
      last_sign=signs[ptr];
      num_changes++;
    }
  return(num_changes>1);
}
/*
int rig_num_zeros_D(int st, int en, int_double *re_zs1, int_double *re_zs2, FILE *outfile, q_state &qs)
{
  int count=0,ptr=st,offset=1,unk_ptr,unk_offset;
  sign_t last_sign,this_sign;
  bool unk_p=false;
  sign_t new_zeros[5];
  int new_zero_ptr=1;

  last_sign=sign(re_zs1[ptr]);
  new_zeros[0]=last_sign;

  if((last_sign==UNK)&&(st==0))
    {
      printf("F(0,chi) is of unknown sign.\n");
      print_int_double_str("F(0,chi)=",re_zs1[0]);
      print_int_double_str("F(1/2B,chi)=",re_zs1[1]);
      print_int_double_str("F(0,chi_bar)=",re_zs2[0]);
      print_int_double_str("F(1/2B,chi_bar)=",re_zs2[1]);
      //qs.type=RE_Z_0_PROB;
      //save_it(qs,outfile);
    }

  // skip initial unknowns
  // assumes region contains at least one known.
  while(last_sign==UNK)
    {
      if(offset==UP_SAMPLE_RATE)
	{
	  offset=1;
	  ptr++;
	  last_sign=sign(re_zs1[ptr]);
	}
      else
	last_sign=check_sign(ptr,offset++,re_zs1,re_zs2);
      new_zeros[new_zero_ptr++]=last_sign;
      if(new_zero_ptr==5)
	new_zero_ptr=1;
    }

  while(ptr<en)
    {
      if(offset==UP_SAMPLE_RATE)
	{
	  offset=1;
	  ptr++;
	  this_sign=sign(re_zs1[ptr]);
	}
      else
	this_sign=check_sign(ptr,offset++,re_zs1,re_zs2);

      // don't populate with unknowns
      if(this_sign==UNK)
	new_zeros[new_zero_ptr++]=last_sign;
      else
	new_zeros[new_zero_ptr++]=this_sign;
      if(new_zero_ptr==5)
	{
	  if(new_zero(new_zeros))
	    {
	      printf("New pair of zeros detected at %d upsample rate.\n",UP_SAMPLE_RATE);
	      printf("Detected at ptr=%d,offset=%d-%d.\n",ptr,offset-5,offset-1);
	    }
	  new_zero_ptr=1;
	  new_zeros[0]=new_zeros[4];
	}

      if(this_sign==UNK)
	{
	  if(!unk_p)
	    {
	      unk_p=true;
	      unk_ptr=ptr;
	      unk_offset=offset-1;
	    }
	  continue;
	}
      if(this_sign==last_sign)
	{
	  if(unk_p)
	    {
	      qs.type=CHECK_PROB;
	      qs.n0=unk_ptr;qs.n02=ptr;
	      qs.offset1=unk_offset;qs.offset2=offset-1;
	      printf("Sign strayed into unknown and back near n0=%d-%d offset=%d-%d.",unk_ptr,ptr,unk_offset,offset-1);
	      if(this_sign==POS)
		{
		  printf(" Look for a negative.\n");
		  qs.exp_sign=NEG;
		}
	      else
		{
		  printf(" Look for a positive.\n");
		  qs.exp_sign=POS;
		}
	      save_it(qs,outfile);
	      unk_p=false;
	    }
	  continue;
	}
      unk_p=false;
      count++;
      last_sign=this_sign;
    }
  return(count);
}
*/
// check for changes of sign between st and en inclusive
// doesn't try to locate missed zeros
int rig_num_zeros_first(int st, int en, int_double *re_zs1, int_double *re_zs2)
{
  int count=0,ptr=st,offset=1;
  sign_t last_sign,this_sign;

  last_sign=sign(re_zs1[ptr]);
  if((last_sign==UNK)&&(st==0))
    {
      printf("F(0,chi) is of unknown sign.\n");
      print_int_double_str("F(0,chi)=",re_zs1[0]);
      print_int_double_str("F(1/2B,chi)=",re_zs1[1]);
      print_int_double_str("F(0,chi_bar)=",re_zs2[0]);
      print_int_double_str("F(1/2B,chi_bar)=",re_zs2[1]);
    }
  // skip initial unknowns
  // assumes region contains at least one known.
  while(last_sign==UNK)
    {
      if(offset==UP_SAMPLE_RATE)
	{
	  offset=1;
	  ptr++;
	  last_sign=sign(re_zs1[ptr]);
	}
      else
	last_sign=check_sign(ptr,offset++,re_zs1,re_zs2);
    }

  while(ptr<en)
    {
      if(offset==UP_SAMPLE_RATE)
	{
	  offset=1;
	  ptr++;
	  this_sign=sign(re_zs1[ptr]);
	}
      else
	this_sign=check_sign(ptr,offset++,re_zs1,re_zs2);
      if((this_sign==last_sign)||(this_sign==UNK))
	continue;
      count++;
      last_sign=this_sign;
    }
  return(count);
}


// check for changes of sign between st and en inclusive
// detects zeros that previous UP_SAMPLE_RATE (/4) failed to spot
int rig_num_zeros_C(int st, int en, int_double *re_zs1, int_double *re_zs2)
{
  int low_count=0,hi_count=0,ptr=st,offset=1;
  sign_t last_low_sign;
  sign_t last_hi_sign,this_hi_sign;
  int lst_chg_off=0,lst_chg_ptr=0;

  last_hi_sign=sign(re_zs1[ptr]);
  last_low_sign=last_hi_sign;


  // is F(0) UNK
  if((last_hi_sign==UNK)&&(st==0))
    {
      printf("F(0,chi) is of unknown sign.\n");
      print_int_double_str("F(0,chi)=",re_zs1[0]);
      print_int_double_str("F(1/2B,chi)=",re_zs1[1]);
      print_int_double_str("F(0,chi_bar)=",re_zs2[0]);
      print_int_double_str("F(1/2B,chi_bar)=",re_zs2[1]);
    }

  // skip initial unknowns
  // assumes region contains at least one known.
  while(last_hi_sign==UNK)
    {
      if(offset==UP_SAMPLE_RATE)
	{
	  offset=1;
	  ptr++;
	  last_hi_sign=sign(re_zs1[ptr]);
	}
      else
	last_hi_sign=check_sign(ptr,offset++,re_zs1,re_zs2);
      if((offset&3)==1)
	last_low_sign=last_hi_sign;
    }

  // start search
  // last_sign=POS/NEG
  while(ptr<en)
    {
      if(offset==UP_SAMPLE_RATE)
	{
	  offset=1;
	  ptr++;
	  this_hi_sign=sign(re_zs1[ptr]);
	}
      else
	this_hi_sign=check_sign(ptr,offset++,re_zs1,re_zs2);

      if(this_hi_sign==UNK)
	continue;
      if(this_hi_sign!=last_hi_sign)
	{
	  hi_count++;
	  last_hi_sign=this_hi_sign;
	}
      if((offset&3)==1)
	if(this_hi_sign!=last_low_sign)
	  {
	    low_count++;
	    last_low_sign=this_hi_sign;
	  }
      if(hi_count>=low_count+2)
	{
	  printf("New pair of zeros detected upsampling at x%d in ptr=%d-%d offset=%d-%d.\n",
		 UP_SAMPLE_RATE,lst_chg_ptr,ptr,lst_chg_off,offset-1);
	  low_count=hi_count;
	}
      lst_chg_ptr=ptr;
      lst_chg_off=offset-1;
    }
  return(hi_count);
}
// check for changes of sign between st and en inclusive
// detects zeros that previous UP_SAMPLE_RATE (/4) failed to spot
// detects strings of ++++?????++++ and ----????----
int rig_num_zeros_D(int st, int en, int_double *re_zs1, int_double *re_zs2, FILE *outfile, q_state &qs)
{
  int low_count=0,hi_count=0,ptr=st,offset=1;
  sign_t last_low_sign;
  sign_t last_hi_sign,this_hi_sign;
  int lst_chg_off=0,lst_chg_ptr=0;
  bool unk_p=false;
  int unk_ptr,unk_offset;

  last_hi_sign=sign(re_zs1[ptr]);
  last_low_sign=last_hi_sign;


  // is F(0) UNK
  if((last_hi_sign==UNK)&&(st==0))
    {
      printf("F(0,chi) is of unknown sign.\n");
      print_int_double_str("F(0,chi)=",re_zs1[0]);
      print_int_double_str("F(1/2B,chi)=",re_zs1[1]);
      print_int_double_str("F(0,chi_bar)=",re_zs2[0]);
      print_int_double_str("F(1/2B,chi_bar)=",re_zs2[1]);
    }

  // skip initial unknowns
  // assumes region contains at least one known.
  while(last_hi_sign==UNK)
    {
      if(offset==UP_SAMPLE_RATE)
	{
	  offset=1;
	  ptr++;
	  last_hi_sign=sign(re_zs1[ptr]);
	}
      else
	last_hi_sign=check_sign(ptr,offset++,re_zs1,re_zs2);
      if((offset&3)==1)
	last_low_sign=last_hi_sign;
    }

  // start search
  // last_sign=POS/NEG
  while(ptr<en)
    {
      if(offset==UP_SAMPLE_RATE)
	{
	  offset=1;
	  ptr++;
	  this_hi_sign=sign(re_zs1[ptr]);
	}
      else
	this_hi_sign=check_sign(ptr,offset++,re_zs1,re_zs2);
      if(this_hi_sign==UNK)
	{
	  if(!unk_p)
	    {
	      unk_p=true;
	      unk_ptr=ptr;
	      unk_offset=offset-1;
	    }
	  continue;
	}
      if(this_hi_sign==last_hi_sign)
	{
	  if(unk_p)
	    {
	      qs.type=CHECK_PROB;
	      qs.n0=unk_ptr;qs.n02=ptr;
	      qs.offset1=unk_offset;qs.offset2=offset-1;
	      printf("Sign strayed into unknown and back near n0=%d-%d offset=%d-%d.",unk_ptr,ptr,unk_offset,offset-1);
	      if(this_hi_sign==POS)
		{
		  printf(" Look for a negative.\n");
		  qs.exp_sign=NEG;
		}
	      else
		{
		  printf(" Look for a positive.\n");
		  qs.exp_sign=POS;
		}
	      save_it(qs,outfile);
	      unk_p=false;
	    }
	}
      else
	{
	  unk_p=false;
	  hi_count++;
	  last_hi_sign=this_hi_sign;
	}

      if((offset&3)==1)
	if(this_hi_sign!=last_low_sign)
	  {
	    low_count++;
	    last_low_sign=this_hi_sign;
	  }

      if(hi_count>=low_count+2)
	{
	  printf("New pair of zeros detected upsampling at x%d in ptr=%d-%d offset=%d-%d.\n",
		 UP_SAMPLE_RATE,lst_chg_ptr,ptr,lst_chg_off,offset-1);
	  low_count=hi_count;
	}

      lst_chg_ptr=ptr;
      lst_chg_off=offset-1;
    }
  return(hi_count);
}

// return the number of zeros in out_array between 0 and UP_SAMPLE_RATE*width inclusive
// the signs of the end points are known to be st_sign and en_sign, neither of which is UNK
// re_z_ptr points at the Z value at the start of the sample
inline int up_num_zeros (mcomplex *out_array, im_s *im_s_vec, sign_t st_sign, sign_t en_sign, int width,
		  int re_z_ptr, int_double *re_zs1, int_double *re_zs2)
{
  int ptr,count=0,last_ptr=0,p;
  sign_t this_sign,last_sign=st_sign;

  for(ptr=1;ptr<=width*UP_SAMPLE_RATE;ptr++)
    {
      this_sign=d_sign(real(out_array[ptr]));
      if((this_sign==last_sign)||(this_sign==UNK))
	continue;
      count++;
      if(count>1)
	break;
      last_sign=this_sign;
    }
  if(count<=1)
    return(count);
  count=0;
  last_sign=st_sign;
  if(width>1)
    {
      for(ptr=1,p=1;ptr<width*UP_SAMPLE_RATE;ptr++,p++)
	{
	  if(p==UP_SAMPLE_RATE)
	    {
	      p=0;
	      re_z_ptr++;
	      this_sign=UNK;
	    }
	  else
	    this_sign=check_sign(re_z_ptr,p,re_zs1,re_zs2);
	  if((this_sign==last_sign)||(this_sign==UNK))
	    continue;
	  count++;
	  last_sign=this_sign;
	}
    }
  else
    {
      for(ptr=1;ptr<width*UP_SAMPLE_RATE;ptr++)
	{
	  this_sign=check_sign(re_z_ptr,ptr,re_zs1,re_zs2);
	  if((this_sign==last_sign)||(this_sign==UNK))
	    continue;
	  count++;
	  last_sign=this_sign;
	}
    }
  if(last_sign!=en_sign)
    count++;
  return(count);
}

inline int up_num_zeros1 (mcomplex *out_array, im_s *im_s_vec, sign_t st_sign, sign_t en_sign, int width,
		   int re_z_ptr, int_double *re_zs1, int_double *re_zs2, FILE *outfile, q_state &qs)
{
  return(0);
}
/*
  int ptr,count=0,last_ptr=0,p,pp;
  sign_t this_sign,last_sign=st_sign;
  int_double this_z;
  for(ptr=1;ptr<=width*UP_SAMPLE_RATE;ptr++)
    {
      this_sign=d_sign(real(out_array[ptr]));
      if((this_sign==last_sign)||(this_sign==UNK))
	continue;
      if(last_ptr==0) // first sign change, this is a given
	{
	  last_sign=this_sign;
	  count++;
	  last_ptr=ptr;
	  continue;
	}
      for(p=last_ptr;p<ptr;p++)
	if(check_sign(re_z_ptr,p,re_zs1,re_zs2)==last_sign)
	  {
	    count++;
	    //last_sign=this_sign;
	    //printf("breaking with p=%d ptr=%d\n",p,ptr);
	    break;
	  }
      if(p==ptr) // failed to find alleged change of sign
	{
	  printf("Check_sign failed to confirm change in sign at n0=%d offset=%d-%d\n",
		 re_z_ptr,last_ptr,ptr-1);
//	  for(int oset=last_ptr;oset<ptr;oset++)
//	    print_int_double_str("Rigorous interval was ",check_sign1(re_z_ptr,oset,re_zs1,re_zs2));
	  qs.type=CHECK_PROB;qs.n0=re_z_ptr;qs.offset1=last_ptr;qs.offset2=ptr-1;qs.exp_sign=last_sign;
	  save_it(qs,outfile);
	}
      last_sign=this_sign;
      last_ptr=ptr;
    }
  return(count);
}
*/

inline void up_sample_load (int_double *re_zs)
{
	for(int i=0;i<UP_SAMPLE_IN_WIDTH-1;i++)
		in_array[i]=avg_int(re_zs[i]);
	up_sample();
}

inline void up_sample_load_neg (int_double *re_zs1, int_double *re_zs2, int ptr)
{
  int i;
  
    if(sign(re_zs1[0])!=sign(re_zs2[0]))
      {
	printf("Error: Signs of re_zs1[0] and re_zs2[0] disgaree.\n");
	print_int_double_str("re_zs1[0]=",re_zs1[0]);
	print_int_double_str("re_zs2[0]=",re_zs2[0]);
      }
/*
    for(i=0;i<UP_SAMPLE_SACRIFICE;i++)
    in_array[i]=-avg_int(re_zs2[UP_SAMPLE_SACRIFICE-i])*exppittwos[UP_SAMPLE_SACRIFICE-i-1];
    else
  */
  for(i=0;i<UP_SAMPLE_SACRIFICE;i++)
    in_array[i]=avg_int(re_zs2[UP_SAMPLE_SACRIFICE-i])*exppittwos[UP_SAMPLE_SACRIFICE-i-1];

  for(;i<UP_SAMPLE_IN_WIDTH-1;i++)
    in_array[i]=avg_int(re_zs1[i-UP_SAMPLE_SACRIFICE]);
  up_sample();
}


// work out the integral of N(T) from re_zs[a] to [b]
// assumes the section a to b fits into the usable
// portion of a single up_sample if re_zs[num_s] goes
// in the end
int_double up_n_twiddle( int a,  int b, int_double *re_zs, im_s *im_s_vec, int t_sam_start)
{
  int i,j,n,out_ptr=t_sam_start,n1=0;
  int_double end_point=int_double(b*one_over_two_B/*im_s_vec[b].im_s*/),res=int_double(0.0);
  sign_t st_sign=sign(re_zs[a]),en_sign;

  while(st_sign==UNK)
    st_sign=sign(re_zs[++a]);
  while(sign(re_zs[b])==UNK)
    b--;
  //printf("Loading from re_zs[%d]\n",out_ptr);
  up_sample_load(&re_zs[out_ptr]);
  //for(i=0;i<UP_SAMPLE_OUT_WIDTH;i++)
  //  printf("%10.8e\n",real(out_array[i]));

  out_ptr=(a-out_ptr)*UP_SAMPLE_RATE;
  //d_print_signs((b-a)*UP_SAMPLE_RATE,&out_array[out_ptr],1);
  //printf("out_ptr=%d\n",out_ptr);
  i=a;
  j=a+1;
  while(j<=b)
    {
      en_sign=sign(re_zs[j]);
      while(en_sign==UNK)
	en_sign=sign(re_zs[++j]);
      n=up_num_zeros(&out_array[out_ptr],
		     im_s_vec,st_sign,en_sign,(j-i),i,re_zs,re_zs); // last arg is dummy here
      if(n!=0)
	{
	  //printf("Found %d zeros here\n",n);
	  //d_print_signs((j-i)*UP_SAMPLE_RATE,&out_array[out_ptr],1);
	  res+=(end_point-int_double(i*one_over_two_B,j*one_over_two_B/*im_s_vec[i].im_s,im_s_vec[j].im_s*/))*n;
	  n1+=n;
	}
      out_ptr+=UP_SAMPLE_RATE*(j-i);
      i=j;
      j++;
      st_sign=en_sign;
    }
  //exit(0);
  //printf("Found %d zeros in Turing Zone.\n",n1);
  return(res);
}

int_double up_n_twiddle1( int a,  int b, int_double *re_zs, im_s *im_s_vec, int t_sam_start, FILE *outfile, q_state &qs)
{
	int i,j,n,out_ptr=t_sam_start,n1=0;
	int_double end_point=int_double(b*one_over_two_B/*im_s_vec[b].im_s*/),res=int_double(0.0);
	sign_t st_sign=sign(re_zs[a]),en_sign;

	while(st_sign==UNK)
		st_sign=sign(re_zs[++a]);
	while(sign(re_zs[b])==UNK)
		b--;
	//printf("Loading from re_zs[%d]\n",out_ptr);
	up_sample_load(&re_zs[out_ptr]);
	//for(i=0;i<UP_SAMPLE_OUT_WIDTH;i++)
	//  printf("%10.8e\n",real(out_array[i]));

	out_ptr=(a-out_ptr)*UP_SAMPLE_RATE;
	//d_print_signs((b-a)*UP_SAMPLE_RATE,&out_array[out_ptr],1);
	//printf("out_ptr=%d\n",out_ptr);
	i=a;
	j=a+1;
	while(j<=b)
	{
		en_sign=sign(re_zs[j]);
		while(en_sign==UNK)
			en_sign=sign(re_zs[++j]);
		n=up_num_zeros1(&out_array[out_ptr],
				im_s_vec,st_sign,en_sign,(j-i),i,re_zs,re_zs,outfile,qs); // last but one arg is dummy here
		if(n!=0)
		{
			//printf("Found %d zeros here\n",n);
			//d_print_signs((j-i)*UP_SAMPLE_RATE,&out_array[out_ptr],1);
		  res+=(end_point-int_double(i*one_over_two_B,j*one_over_two_B/*im_s_vec[i].im_s,im_s_vec[j].im_s*/))*n;
			n1+=n;
		}
		out_ptr+=UP_SAMPLE_RATE*(j-i);
		i=j;
		j++;
		st_sign=en_sign;
	}
	//exit(0);
	//printf("Found %d zeros in Turing Zone.\n",n1);
	return(res);
}

// originally used Rumely's Th2 p 429
// int S(t) <= 1.8397+0.1242*log(q*t2/2/Pi)
//
// now use Trudgian's optimised for qt=10^8
//          <= 2.17618+0.0679955*log(q*t2/2/Pi)
//
inline int_double s_of_t( int q, double t2)
{
  int_double tmp2=int_double (t2/2.0),tmp=int_double(67996), res=int_double(217618);
  tmp/=1000000;res/=100000;
  tmp2*=q;
  tmp2/=d_pi;
  res=res+tmp*log(tmp2);
  res.left=res.right;
  return(res);
}

// returns (t1+t2)*ln(q/pi)/4
inline int_double ln_term (double t1, double t2,  int q)
{
	return(log(int_double(q)/d_pi)*(t1+t2)/4);
}


sign_t up_sam(int ptr, int_double delta, int_double *re_zs)
{
	int_double res=d_inter_err;
	int_double del;
	int n;
	print_int_double_str("up_sam called with delta =",delta);
	printf("ptr= %d.\n",ptr);
	del=delta-one_over_two_B;
	for(n=1;n<=INTER_N;n++)
	{
		res+=re_zs[ptr+n]*exp(-del*del/TWO_H_SQR)*sincpi(twoBpi*del);
		print_int_double_str("re_zs= ",re_zs[ptr+n]);
		print_int_double_str("res  = ",res);
		del-=one_over_two_B;
	}
	del=delta;
	for(n=0;n<INTER_N;n++)
	{
		res+=re_zs[ptr-n]*exp(-del*del/TWO_H_SQR)*sincpi(twoBpi*del);
		print_int_double_str("re_zs= ",re_zs[ptr+n]);
		print_int_double_str("res  = ",res);
		del+=one_over_two_B;
	}
	print_int_double_str("up_sam computed ",res);
	return(sign(res));
}
/*
inline bool check_sign1(int ptr, int_double x, sign_t exp_sign, int_double *re_zs)
{
//printf("In check_sign1 with ptr=%d\n.",ptr);
if(x.right>0.0) // negative
{
ptr--;
x+=1;
}
// now our point is between re_zs[ptr] and re_zs[ptr+1]
//return(up_sam(ptr,x*one_over_two_B,re_zs)==exp_sign);
return(check_sign(ptr,x.left*UP_SAMPLE_RATE,re_zs,re_zs)==exp_sign);
}
*/

// count zeros between a and b
inline int num_zeros1(int a, int b, int_double *re_zs)
{
	int zero_count=0,ptr=a+1;
	sign_t this_sign,last_sign=sign(re_zs[a]);
	while(last_sign==UNK)
	{
		last_sign=sign(re_zs[ptr++]);
		if(ptr>b)
			fatal_error("Catastrophic failure in num_zeros. Exiting.\n");
	}

	while(ptr<=b)
	{
		this_sign=sign(re_zs[ptr++]);
		while((this_sign==UNK)&&(ptr<=b))
			this_sign=sign(re_zs[ptr++]);
		if(this_sign==UNK) // file finished with one or more unknowns
			return(zero_count);
		if(this_sign!=last_sign)
		{
			zero_count++;
			last_sign=this_sign;
		}
	}
	return(zero_count);
}	

int stat_count=0;

inline int_double calc_int (double st, double en, double endpt)
{
  return(int_double(endpt-en,endpt-st));
}

int_double rig_n_twiddle1 (int t1, int t2, int_double *re_zs, im_s *im_s_vec, FILE *outfile, q_state &qs)
{
  double delta=(one_over_two_B/*(im_s_vec[t2].im_s-im_s_vec[t2-1].im_s*/)/(double) UP_SAMPLE_RATE,endpt=t2*one_over_two_B;//im_s_vec[t2].im_s;
  double last_change=t1*one_over_two_B/*im_s_vec[t1].im_s*/,this_change;
  int_double res=d_zero;
  int ptr=t1,offset=1;
  sign_t last_sign=sign(re_zs[t1]),this_sign;
  bool unk_p=false; int unk_ptr,unk_offset;
    while(ptr<t2)
    {
      if(offset==UP_SAMPLE_RATE)
	{
	  offset=1;
	  ptr++;
	  //printf("ptr=%d\n",ptr);
	  this_sign=sign(re_zs[ptr]);
	}
      else
	this_sign=check_sign(ptr,offset++,re_zs,re_zs);
      if(this_sign==last_sign)
	{
	  if(unk_p)
	    {
	      qs.type=CHECK_PROB;
	      qs.n0=unk_ptr;qs.n02=ptr;
	      qs.offset1=unk_offset;qs.offset2=offset-1;
	      printf("Sign strayed into unknown and back in TZ near n0=%d-%d offset=%d-%d.",unk_ptr,ptr,unk_offset,offset-1);
	      if(this_sign==POS)
		{
		  printf(" Look for a negative.\n");
		  qs.exp_sign=NEG;
		}
	      else
		{
		  printf(" Look for a positive.\n");
		  qs.exp_sign=POS;
		}
	      save_it(qs,outfile);
	      unk_p=false;
	    }
	  last_change+=delta;
	  continue;
	}
      if(this_sign==UNK)
	{
	  if(!unk_p)
	    {
	      unk_p=true;
	      unk_ptr=ptr;
	      unk_offset=offset-1;
	    }
	  continue;
	}
      unk_p=false;
      this_change=ptr*one_over_two_B/*im_s_vec[ptr].im_s*/+(offset-1)*delta;
      //printf("Zero found at ptr=%d between offset=%d and offset=%d\n",ptr,offset-2,offset-1);
      res+=calc_int(last_change,this_change,endpt);
      last_change=this_change;
      //print_int_double_str("res is now=",res);
      last_sign=this_sign;
    }
    return(res);
}

int_double rig_n_twiddle (int t1, int t2, int_double *re_zs, im_s *im_s_vec)
{
  double delta=(one_over_two_B/*im_s_vec[t2].im_s-im_s_vec[t2-1].im_s*/)/(double) UP_SAMPLE_RATE,endpt=t2*one_over_two_B;//im_s_vec[t2].im_s;
  double last_change=t1*one_over_two_B/*im_s_vec[t1].im_s*/,this_change;
  int_double res=d_zero;
  int ptr=t1,offset=1;
  sign_t last_sign=sign(re_zs[t1]),this_sign;
  bool unk_p=false; int unk_ptr,unk_offset;
    while(ptr<t2)
    {
      if(offset==UP_SAMPLE_RATE)
	{
	  offset=1;
	  ptr++;
	  //printf("ptr=%d\n",ptr);
	  this_sign=sign(re_zs[ptr]);
	}
      else
	this_sign=check_sign(ptr,offset++,re_zs,re_zs);
      if(this_sign==last_sign)
	{
	  if(unk_p)
	    {
	      //qs.type=CHECK_PROB;
	      //qs.n0=unk_ptr;qs.n02=ptr;
	      //qs.offset1=unk_offest;qs.offset2=offset-1;
	      printf("Sign strayed into unknown and back in TZ near n0=%d-%d offset=%d-%d.",unk_ptr,ptr,unk_offset,offset-1);
	      if(this_sign==POS)
		{
		  printf(" Look for a negative.\n");
		  //qs.exp_sign=NEG;
		}
	      else
		{
		  printf(" Look for a positive.\n");
		  //qs.exp_sign=POS;
		}
	      //save_it(qs,outfile);
	      unk_p=false;
	    }
	  last_change+=delta;
	  continue;
	}
      if(this_sign==UNK)
	{
	  if(!unk_p)
	    {
	      unk_p=true;
	      unk_ptr=ptr;
	      unk_offset=offset-1;
	    }
	  continue;
	}
      unk_p=false;
      this_change=ptr*one_over_two_B/*im_s_vec[ptr].im_s*/+(offset-1)*delta;
      //printf("Zero found at ptr=%d between offset=%d and offset=%d\n",ptr,offset-2,offset-1);
      res+=calc_int(last_change,this_change,endpt);
      last_change=this_change;
      //print_int_double_str("res is now=",res);
      last_sign=this_sign;
    }
    return(res);
}

void qs_omega(q_state &qs, int_complex &omega2)
{
  qs.omega[0]=omega2.real.left;
  qs.omega[1]=-omega2.real.right;
  qs.omega[2]=omega2.imag.left;
  qs.omega[3]=-omega2.imag.right;
  qs.index=-qs.index;
}

int rig_calc_zeros_cmplx1( int t1_ptr, int t2_ptr, int q,
			  int_double &gamma_int, int_double *re_zs1, int_double *re_zs2, im_s *im_s_vec,
			   int_complex &omega2,
			   int index, FILE *outfile, q_state &qs)
{
  //printf("Reached rig_calc_zeros.\n");
  //exit(0);
  double t1=t1_ptr*one_over_two_B;//im_s_vec[t1_ptr].im_s;
  double t2=t2_ptr*one_over_two_B;//im_s_vec[t2_ptr].im_s;
  int_double ntwiddle1=rig_n_twiddle1(t1_ptr,t2_ptr,re_zs1,im_s_vec,outfile,qs);
  qs_omega(qs,omega2);
  int_double ntwiddle2=rig_n_twiddle1(t1_ptr,t2_ptr,re_zs2,im_s_vec,outfile,qs);
  
  int_double s_t=s_of_t(q,t2)*2;
	int_double lnt=ln_term(t1,t2,q)*2;
	int_double calc_zeros;
	int calc_zeros_int;

	calc_zeros=(lnt+gamma_int*2/(t2-t1))/d_pi-(ntwiddle1+ntwiddle2)/(t2-t1)+s_t/(t2-t1); // Turing's estimate
	//print_int_double_str("Turing's Estimate = ",calc_zeros);

	calc_zeros_int=(int)ceil(calc_zeros.left);
	//printf("Turing's estimate = %d\n",calc_zeros_int);

	if(calc_zeros_int!=floor(-calc_zeros.right)) // Turing's estimate did not bracket a unique integer
	{
		printf("q: %ld index: %ld Problem with Turing's estimate:- ",q,index);
		print_int_double(calc_zeros);
		printf("\n");
		return(BAD_DIFF);
	}
	return(calc_zeros_int);
}
int rig_calc_zeros_cmplx( int t1_ptr, int t2_ptr, int q,
			  int_double &gamma_int, int_double *re_zs1, int_double *re_zs2, im_s *im_s_vec,
				int index)
{
  //printf("Reached rig_calc_zeros.\n");
  //exit(0);
  double t1=t1_ptr*one_over_two_B;//im_s_vec[t1_ptr].im_s;
  double t2=t2_ptr*one_over_two_B;//im_s_vec[t2_ptr].im_s;
  int_double ntwiddle1=rig_n_twiddle(t1_ptr,t2_ptr,re_zs1,im_s_vec);
  int_double ntwiddle2=rig_n_twiddle(t1_ptr,t2_ptr,re_zs2,im_s_vec);
  print_int_double_str("n1=",ntwiddle1);
  print_int_double_str("n2=",ntwiddle2);

  int_double s_t=s_of_t(q,t2)*2;
  int_double lnt=ln_term(t1,t2,q)*2;
  int_double calc_zeros;
  int calc_zeros_int;
  print_int_double_str("lnt=",lnt);
  print_int_double_str("s(t)=",s_t);

  calc_zeros=(lnt+gamma_int*2/(t2-t1))/d_pi-(ntwiddle1+ntwiddle2)/(t2-t1)+s_t/(t2-t1); // Turing's estimate
	//print_int_double_str("Turing's Estimate = ",calc_zeros);

  calc_zeros_int=(int)ceil(calc_zeros.left);
	//printf("Turing's estimate = %d\n",calc_zeros_int);

  if(calc_zeros_int!=floor(-calc_zeros.right)) // Turing's estimate did not bracket a unique integer
    {
      printf("q: %ld index: %ld Problem with Turing's estimate:- ",q,index);
      print_int_double(calc_zeros);
      printf("\n");
      return(BAD_DIFF);
    }
  return(calc_zeros_int);
}

int rig_calc_zeros_re1( int t1_ptr, int t2_ptr, int q,
				int_double &gamma_int, int_double &arg_omega,int_double *re_zs, im_s *im_s_vec,
			int index, FILE *outfile, q_state &qs)
{
  //printf("Reached rig_calc_zeros.\n");
  //exit(0);
  double t1=t1_ptr*one_over_two_B;//im_s_vec[t1_ptr].im_s;
  double t2=t2_ptr*one_over_two_B;//im_s_vec[t2_ptr].im_s;
	int_double ntwiddle=rig_n_twiddle1(t1_ptr,t2_ptr,re_zs,im_s_vec,outfile,qs);

	int_double s_t=s_of_t(q,t2);
	int_double lnt=ln_term(t1,t2,q);
	int_double calc_zeros;
	int calc_zeros_int;

	calc_zeros=(arg_omega+lnt+gamma_int/(t2-t1))/d_pi-ntwiddle/(t2-t1)+s_t/(t2-t1); // Turing's estimate
	//print_int_double_str("Turing's Estimate = ",calc_zeros);

	calc_zeros_int=(int)ceil(calc_zeros.left);
	printf("Turing's estimate = %d\n",calc_zeros_int);

	if(calc_zeros_int!=floor(-calc_zeros.right)) // Turing's estimate did not bracket a unique integer
	{
		printf("q: %ld index: %ld Problem with Turing's estimate:- ",q,index);
		print_int_double(calc_zeros);
		printf("\n");
		return(BAD_DIFF);
	}
	return(calc_zeros_int);
}
int rig_calc_zeros_re( int t1_ptr, int t2_ptr, int q,
				int_double &gamma_int, int_double &arg_omega,int_double *re_zs, im_s *im_s_vec,
				int index)
{
  //printf("Reached rig_calc_zeros.\n");
  //exit(0);
  double t1=t1_ptr*one_over_two_B;//im_s_vec[t1_ptr].im_s;
  double t2=t2_ptr*one_over_two_B;//im_s_vec[t2_ptr].im_s;
	int_double ntwiddle=rig_n_twiddle(t1_ptr,t2_ptr,re_zs,im_s_vec);

	int_double s_t=s_of_t(q,t2);
	int_double lnt=ln_term(t1,t2,q);
	int_double calc_zeros;
	int calc_zeros_int;

	calc_zeros=(arg_omega+lnt+gamma_int/(t2-t1))/d_pi-ntwiddle/(t2-t1)+s_t/(t2-t1); // Turing's estimate
	//print_int_double_str("Turing's Estimate = ",calc_zeros);

	calc_zeros_int=(int)ceil(calc_zeros.left);
	//printf("Turing's estimate = %d\n",calc_zeros_int);

	if(calc_zeros_int!=floor(-calc_zeros.right)) // Turing's estimate did not bracket a unique integer
	{
		printf("q: %ld index: %ld Problem with Turing's estimate:- ",q,index);
		print_int_double(calc_zeros);
		printf("\n");
		return(BAD_DIFF);
	}
	return(calc_zeros_int);
}


// compare number of zeros observed with that calculate by Turing's method
int check_zeros(int count_zeros, int t1_ptr, int t2_ptr, int q,
				int_double &gamma_int, int_double &arg_omega,int_double *re_zs, im_s *im_s_vec,
				int index, int num_s)
{
  double t1=t1_ptr*one_over_two_B;//im_s_vec[t1_ptr].im_s;
  double t2=t2_ptr*one_over_two_B;//im_s_vec[t2_ptr].im_s;
	int_double ntwiddle=up_n_twiddle(t1_ptr,t2_ptr,re_zs,im_s_vec,num_s);

	int_double s_t=s_of_t(q,t2);
	int_double lnt=ln_term(t1,t2,q);
	int_double calc_zeros;
	int calc_zeros_int;

	calc_zeros=(arg_omega+lnt+gamma_int/(t2-t1))/d_pi-ntwiddle/(t2-t1)+s_t/(t2-t1); // Turing's estimate
	//print_int_double_str("Turing's Estimate = ",calc_zeros);

	calc_zeros_int=(int)ceil(calc_zeros.left);
	//printf("Turing's estimate = %d\n",calc_zeros_int);

	if(calc_zeros_int!=floor(-calc_zeros.right)) // Turing's estimate did not bracket a unique integer
	{
		printf("q: %ld index: %ld Problem with Turing's estimate:- ",q,index);
		print_int_double(calc_zeros);
		printf("\n");
		return(BAD_DIFF);
	}
	return(calc_zeros_int-count_zeros);
}

int check_zeros1(int count_zeros, int t1_ptr, int t2_ptr, int q,
				int_double &gamma_int, int_double &arg_omega,int_double *re_zs, im_s *im_s_vec,
		 int index, int num_s, FILE *outfile, q_state &qs)
{
  double t1=t1_ptr*one_over_two_B;//im_s_vec[t1_ptr].im_s;
  double t2=t2_ptr*one_over_two_B;//im_s_vec[t2_ptr].im_s;
	int_double ntwiddle=up_n_twiddle1(t1_ptr,t2_ptr,re_zs,im_s_vec,num_s,outfile,qs);

	int_double s_t=s_of_t(q,t2);
	int_double lnt=ln_term(t1,t2,q);
	int_double calc_zeros;
	int calc_zeros_int;

	calc_zeros=(arg_omega+lnt+gamma_int/(t2-t1))/d_pi-ntwiddle/(t2-t1)+s_t/(t2-t1); // Turing's estimate

	calc_zeros_int=(int)ceil(calc_zeros.left);

	if(calc_zeros_int!=floor(-calc_zeros.right)) // Turing's estimate did not bracket a unique integer
	{
		printf("q: %ld index: %ld Problem with Turing's estimate:- ",q,index);
		print_int_double(calc_zeros);
		printf("\n");
		return(BAD_DIFF);
	}
	return(calc_zeros_int-count_zeros);
}

int find_more_zeros(int_double *re_zs,  int_double *re_zs2, int turing_start,  int missing, im_s *im_s_vec)
{
  int n,up_n;
  int t2=turing_start;
  int t1=t2-UP_SAMPLE_USABLE_WIDTH+1;
  int t3=t2+UP_SAMPLE_SACRIFICE;
  int t0=t3-UP_SAMPLE_IN_WIDTH+1;
  //double d1,d2,d3;
  int still_missing=missing;
  int ptr_hi,ptr_lo;
  sign_t st_sign,en_sign;
  while(t0>0)
    {
      up_sample_load(&re_zs[t0]);

      while(sign(re_zs[t1])==UNK) t1--;
      en_sign=sign(re_zs[t2]);
      while(en_sign==UNK) en_sign=sign(re_zs[++t2]);
      ptr_hi=t2;
      ptr_lo=ptr_hi-1;
      while(ptr_hi>t1)
	{
	  st_sign=sign(re_zs[ptr_lo]);
	  while(st_sign==UNK)
	    st_sign=sign(re_zs[--ptr_lo]);
	  up_n=up_num_zeros(&out_array[(ptr_lo-t0)*UP_SAMPLE_RATE],im_s_vec,st_sign,en_sign,
			    ptr_hi-ptr_lo,ptr_lo,re_zs,re_zs2);
	  //if(up_n>1)
	  //  printf("up_num_zeros returned %d\n",up_n);
	  still_missing-=(up_n&0xFFFFFFFE); // if odd, will already know about one of them
	  en_sign=st_sign;
	  ptr_hi=ptr_lo;
	  ptr_lo--;
	}
      t2=t1;
      t1=t2-UP_SAMPLE_USABLE_WIDTH+1;
      t3=t2+UP_SAMPLE_SACRIFICE;
      t0=t3-UP_SAMPLE_IN_WIDTH+1;
    }
  // if here then t0<0 so need to load the negative portion
  // we have checked from t2 upwards so now check
  // from 0 to t2-1
  // we have en_sign=sign of t2
  up_sample_load_neg(re_zs,re_zs2,t2);
  ptr_hi=t2;
  ptr_lo=ptr_hi-1;
  while(ptr_hi>0)
    {
      st_sign=sign(re_zs[ptr_lo]);
      while(st_sign==UNK)
	{
	  ptr_lo--;
	  if(ptr_lo<0)
	    {
	      printf("re_zs[0] is of unknown sign. Cannot resolve it.\n");
	      print_int_double_str("re_zs[0]=",re_zs[0]);
	      print_int_double_str("re_zs[1]=",re_zs[1]);
	      print_int_double_str("re_zs2[0]=",re_zs2[0]);
	      print_int_double_str("re_zs2[1]=",re_zs2[1]);
	      return(still_missing);
	    }
	  st_sign=sign(re_zs[ptr_lo]);
	}
      up_n=up_num_zeros(&out_array[(ptr_lo+UP_SAMPLE_SACRIFICE)*UP_SAMPLE_RATE],im_s_vec,
			st_sign,en_sign,ptr_hi-ptr_lo,ptr_lo,re_zs,re_zs2);
      // if up_n is odd, we will have counted one already
      //if(up_n>1)
      //  printf("up_num_zeros returned %d\n",up_n);

      still_missing-=(up_n&0xFFFFFFFE);

      //if(still_missing==0)
      //return(0);

      en_sign=st_sign;
      ptr_hi=ptr_lo;
      ptr_lo--;
    }
  return(still_missing);
}


int find_more_zeros1(int_double *re_zs,  int_double *re_zs2, int turing_start,  int missing, im_s *im_s_vec, FILE *outfile, q_state &qs)
{
	int n,up_n;
	int t2=turing_start;
	int t1=t2-UP_SAMPLE_USABLE_WIDTH+1;
	int t3=t2+UP_SAMPLE_SACRIFICE;
	int t0=t3-UP_SAMPLE_IN_WIDTH+1;
	//double d1,d2,d3;
	int still_missing=missing;
	int ptr_hi,ptr_lo;
	sign_t st_sign,en_sign;
	while(t0>0)
	{
		up_sample_load(&re_zs[t0]);

		while(sign(re_zs[t1])==UNK) t1--;
		en_sign=sign(re_zs[t2]);
		while(en_sign==UNK) en_sign=sign(re_zs[++t2]);
		ptr_hi=t2;
		ptr_lo=ptr_hi-1;
		while(ptr_hi>t1)
		{
			st_sign=sign(re_zs[ptr_lo]);
			while(st_sign==UNK)
				st_sign=sign(re_zs[--ptr_lo]);
			up_n=up_num_zeros1(&out_array[(ptr_lo-t0)*UP_SAMPLE_RATE],im_s_vec,st_sign,en_sign,
				ptr_hi-ptr_lo,ptr_lo,re_zs,re_zs2,outfile,qs);
			still_missing-=(up_n&0xFFFFFFFE); // if odd, will already know about one of them
			en_sign=st_sign;
			ptr_hi=ptr_lo;
			ptr_lo--;
		}
		t2=t1;
		t1=t2-UP_SAMPLE_USABLE_WIDTH+1;
		t3=t2+UP_SAMPLE_SACRIFICE;
		t0=t3-UP_SAMPLE_IN_WIDTH+1;
	}
	// if here then t0<0 so need to load the negative portion
	// we have checked from t2 upwards so now check
	// from 0 to t2-1
	// we have en_sign=sign of t2
	up_sample_load_neg(re_zs,re_zs2,t2);
	ptr_hi=t2;
	ptr_lo=ptr_hi-1;
	while(ptr_hi>0)
	{
		st_sign=sign(re_zs[ptr_lo]);
		while(st_sign==UNK)
		{
			ptr_lo--;
			if(ptr_lo<0)
			{
				printf("re_zs[0] is of unknown sign. Cannot resolve it.\n");
				print_int_double_str("re_zs[0]=",re_zs[0]);
				print_int_double_str("re_zs[1]=",re_zs[1]);
				print_int_double_str("re_zs2[0]=",re_zs2[0]);
				print_int_double_str("re_zs2[1]=",re_zs2[1]);
				return(still_missing);
			}
			st_sign=sign(re_zs[ptr_lo]);
		}
		up_n=up_num_zeros1(&out_array[(ptr_lo+UP_SAMPLE_SACRIFICE)*UP_SAMPLE_RATE],im_s_vec,
			st_sign,en_sign,ptr_hi-ptr_lo,ptr_lo,re_zs,re_zs2,outfile,qs);
		//if(up_n>1)
		//printf("up_n=%d\n",up_n);
		// if up_n is odd, we will have counted one already
		still_missing-=(up_n&0xFFFFFFFE);

		//if(still_missing==0)
		//return(0);

		en_sign=st_sign;
		ptr_hi=ptr_lo;
		ptr_lo--;
	}
	return(still_missing);
}

inline int find_more_zeros2(int_double *re_zs1, int_double *re_zs2,  int turing_start,  int missing,
							im_s *im_s_vec)
{
	int still_missing=missing;
	still_missing=find_more_zeros(re_zs1,re_zs2,turing_start,missing,im_s_vec);
	if(still_missing!=0)
		still_missing=find_more_zeros(re_zs2,re_zs1,turing_start,still_missing,im_s_vec);
	return(still_missing);
}

bool save_p=false;

void save_state(int q, int index, bool real_p, bool neg_one, int_complex &omega, 
				bool neg_one2, int_complex & omega2, int num_zeros, int num_s, 
				int_double *re_zs1, int_double *re_zs2, FILE *out_file)
{
	save_p=true;
	fwrite(&q,sizeof(int),1,out_file);
	fwrite(&index,sizeof(int),1,out_file);
	fwrite(&real_p,sizeof(bool),1,out_file);
	fwrite(&num_zeros,sizeof(int),1,out_file);
	fwrite(&neg_one,sizeof(bool),1,out_file);
	fwrite(&omega,sizeof(int_complex),1,out_file);
	//fwrite(&num_s,sizeof(int),1,out_file);
	fwrite(re_zs1,sizeof(int_double),num_s,out_file);
	if(!real_p)
	{
		fwrite(&neg_one2,sizeof(bool),1,out_file);
		fwrite(&omega2,sizeof(int_complex),1,out_file);
		fwrite(re_zs2,sizeof(int_double),num_s,out_file);
	}
}

void save_state1(unsigned int q,
		 unsigned int index,
		 bool real_p,
		 bool neg_one,
		 const int_complex &omega,
		 const int_complex &omega2,
		 unsigned int N0,
		 int_double *re_zs1,
		 int_double *re_zs2,
		 unsigned int t_sam_start,
		 unsigned int t_sam_end,
		 const int_double &gamma,
		 FILE *outfile)
{
  fwrite(&q,sizeof(unsigned int),1,outfile);
  fwrite(&index,sizeof(unsigned int),1,outfile);
  fwrite(&real_p,sizeof(bool),1,outfile);
  fwrite(&neg_one,sizeof(bool),1,outfile);
  fwrite(&omega,sizeof(int_complex),1,outfile);
  fwrite(&N0,sizeof(unsigned int),1,outfile);
  fwrite(re_zs1,sizeof(int_double),N0,outfile);
  fwrite(&t_sam_start,sizeof(unsigned int),1,outfile);
  fwrite(&t_sam_end,sizeof(unsigned int),1,outfile);
  fwrite(&gamma,sizeof(int_double),1,outfile);
  if(!real_p)
    {
      fwrite(&omega2,sizeof(int_complex),1,outfile);
      fwrite(re_zs2,sizeof(int_double),N0,outfile);
    }
}


void read_check( int expect,  int i, FILE **infiles)
{
	int tmp;
	if(fread(&tmp,sizeof(int),1,infiles[i]))
		if(tmp==expect)
			return;
	fatal_error("Data mismatch between input files. Exiting.\n");
}


inline int conj_j ( int j,  int num_chi,  int q)
{
	if(q&7) // q not divisible by 8
		return(num_chi-j-1);
	if(j<(num_chi>>1))
		return((num_chi>>1)-j-1);
	return(num_chi+(num_chi>>1)-j-1);
}

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
