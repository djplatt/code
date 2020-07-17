#ifndef _GDD_CU_
#define _GDD_CU_


//#include "gqd_api.h"
#include "gdd_inline.cu"

/** constants in the constant memory */
#define n_dd_inv_fact (15)
static __device__ __constant__ GPU_dd dd_inv_fact[n_dd_inv_fact];
GPU_dd* d_dd_sin_table = NULL;
GPU_dd* d_dd_cos_table = NULL;

/** init function */
void GDDFree() {
	if(d_dd_sin_table) {
		GPUFREE(d_dd_sin_table);
	}

	if(d_dd_cos_table) {
		GPUFREE(d_dd_cos_table);
	}
}

void GDDInit() {
	printf("GDD initialization...\n");
	
	//inverse table
	GPU_dd h_inv_fact[] = {
	  make_dd( 1.66666666666666657e-01,  9.25185853854297066e-18),
	  make_dd( 4.16666666666666644e-02,  2.31296463463574266e-18),
	  make_dd( 8.33333333333333322e-03,  1.15648231731787138e-19),
	  make_dd( 1.38888888888888894e-03, -5.30054395437357706e-20),
	  make_dd( 1.98412698412698413e-04,  1.72095582934207053e-22),
	  make_dd( 2.48015873015873016e-05,  2.15119478667758816e-23),
	  make_dd( 2.75573192239858925e-06, -1.85839327404647208e-22),
	  make_dd( 2.75573192239858883e-07,  2.37677146222502973e-23),
	  make_dd( 2.50521083854417202e-08, -1.44881407093591197e-24),
	  make_dd( 2.08767569878681002e-09, -1.20734505911325997e-25),
	  make_dd( 1.60590438368216133e-10,  1.25852945887520981e-26),
	  make_dd( 1.14707455977297245e-11,  2.06555127528307454e-28),
	  make_dd( 7.64716373181981641e-13,  7.03872877733453001e-30),
	  make_dd( 4.77947733238738525e-14,  4.39920548583408126e-31),
	  make_dd( 2.81145725434552060e-15,  1.65088427308614326e-31)

	};
	CUDA_SAFE_CALL( cudaMemcpyToSymbol( dd_inv_fact, h_inv_fact, sizeof(GPU_dd)*n_dd_inv_fact ) );

	GPU_dd h_sin_table [] = {
	  make_dd(1.950903220161282758e-01, -7.991079068461731263e-18), 
	  make_dd(3.826834323650897818e-01, -1.005077269646158761e-17), 
	  make_dd(5.555702330196021776e-01,  4.709410940561676821e-17),
	  make_dd(7.071067811865475727e-01, -4.833646656726456726e-17)
	};
	GPUMALLOC((void**)&d_dd_sin_table, sizeof(GPU_dd)*4);
	TOGPU(d_dd_sin_table, h_sin_table, sizeof(GPU_dd)*4);

	GPU_dd h_cos_table [] = {
	  make_dd(9.807852804032304306e-01, 1.854693999782500573e-17),
	  make_dd(9.238795325112867385e-01, 1.764504708433667706e-17),
	  make_dd(8.314696123025452357e-01, 1.407385698472802389e-18),
	  make_dd(7.071067811865475727e-01, -4.833646656726456726e-17)
	};
	GPUMALLOC((void**)&d_dd_cos_table, sizeof(GPU_dd)*4);
	TOGPU(d_dd_cos_table, h_cos_table, sizeof(GPU_dd)*4);
}



__device__
GPU_dd exp(const GPU_dd &a) 
{	

	const double k = 512.0;
	const double inv_k = 1.0 / k;

	if (a.x <= -709.0)
		return make_dd(0.0);

	//!!!!!!!!!!!!!
	if (a.x >=  709.0)
		return make_dd(0.0);
	//return dd_real::_inf;

	if (is_zero(a))
		return make_dd(1.0);

	if (is_one(a))
		return _dd_e;

	double m = floor(a.x / _dd_log2.x + 0.5);
	GPU_dd r = mul_pwr2(a - _dd_log2 * m, inv_k);
	GPU_dd s, t, p;

	p = sqr(r);
	s = r + mul_pwr2(p, 0.5);
	p = p * r;
	t = p * dd_inv_fact[0];
	int i = 0;
	do {
		s = s + t;
		p = p * r;
		t = p * dd_inv_fact[++i];
	} while ((fabs(to_double(t)) > inv_k * _dd_eps) && (i < 5));

	s = s + t;

	for( int i = 0; i < 9; i++ )
	{
		s = mul_pwr2(s, 2.0) + sqr(s);
	}

	s = s + 1.0;

//#ifdef NATIVE_DOUBLE
//	return ldexp(s, __double2int_rn(m));
//#else
	return ldexp(s, int(m));
//#endif
}


/* Computes the square root of the double-double number dd.
   NOTE: dd must be a non-negative number.                   */
__device__
GPU_dd sqrt(const GPU_dd &a) 
{
  if (is_zero(a))
    return make_dd(0.0);

  //!!!!!!!!!!!!!!
  //TO DO: should make an error
  if (is_negative(a)) {
    //return _nan;
	  return make_dd( 0.0 );
  }

  double x = 1.0 / sqrt(a.x);
  double ax = a.x * x;
 
  return dd_add(ax, (a - sqr(ax)).x * (x * 0.5));
  //return a - sqr(ax);
}

/* Logarithm.  Computes log(x) in double-double precision.
   This is a natural logarithm (i.e., base e).            */
__device__
GPU_dd log(const GPU_dd &a) 
{
  if (is_one(a)) {
    return make_dd(0.0);
  }

  //!!!!!!!!!
  //TO DO: return an errro
  if (a.x <= 0.0) 
  {
    //return _nan;
	  return make_dd( 0.0 );
  }

  GPU_dd x = make_dd(log(a.x));   // Initial approximation 

  x = x + a * exp(negative(x)) - 1.0;

  return x;
}


/* Computes sin(a) using Taylor series.
   Assumes |a| <= pi/32.                           */
__device__
GPU_dd sin_taylor(const GPU_dd &a) {
  const double thresh = 0.5 * fabs(to_double(a)) * _dd_eps;
  GPU_dd r, s, t, x;

  if (is_zero(a)) {
    return make_dd(0.0);
  }

  int i = 0;
  x = negative(sqr(a)); //-sqr(a);
  s = a;
  r = a;
  do {
    r = r*x;
    t = r * dd_inv_fact[i];
    s = s + t;
    i += 2;
  } while (i < n_dd_inv_fact && fabs(to_double(t)) > thresh);

  return s;
}

__device__
GPU_dd cos_taylor(const GPU_dd &a) {
  const double thresh = 0.5 * _dd_eps;
  GPU_dd r, s, t, x;

  if (is_zero(a)) {
    return make_dd(1.0);
  }

  x = negative(sqr(a));
  r = x;
  s = 1.0 + mul_pwr2(r, 0.5);
  int i = 1;
  do {
    r = r*x;
    t = r * dd_inv_fact[i];
    s = s + t;
    i += 2;
  } while (i < n_dd_inv_fact && fabs(to_double(t)) > thresh);

  return s;
}

__device__
void sincos_taylor(const GPU_dd &a, 
                          GPU_dd &sin_a, GPU_dd &cos_a) {
  if (is_zero(a)) {
    sin_a.x = 0.0; sin_a.y = 0.0;
    cos_a.x = 1.0; cos_a.y = 0.0;
    return;
  }

  sin_a = sin_taylor(a);
  cos_a = sqrt(1.0 - sqr(sin_a));
}


__device__
GPU_dd sin(const GPU_dd &a, const GPU_dd* d_dd_sin_table, const GPU_dd* d_dd_cos_table) {  

  if (is_zero(a)) {
    return make_dd(0.0);
  }

  // approximately reduce modulo 2*pi
  GPU_dd z = nint(a / _dd_2pi);
  GPU_dd r = a - _dd_2pi * z;

  // approximately reduce modulo pi/2 and then modulo pi/16.
  GPU_dd t;
  double q = floor(r.x / _dd_pi2.x + 0.5);
  t = r - _dd_pi2 * q;
  int j = (int)(q);
  q = floor(t.x / _dd_pi16.x + 0.5);
  t = t - _dd_pi16 * q;
  int k = (int)(q);
  int abs_k = abs(k);

  if (j < -2 || j > 2) {
    //dd_real::error("(dd_real::sin): Cannot reduce modulo pi/2.");
    r.x = r.y = 0.0;
    return r;
  }

  if (abs_k > 4) {
    //dd_real::error("(dd_real::sin): Cannot reduce modulo pi/16.");
    r.x = r.y = 0.0;
    return r;
  }

  if (k == 0) {
    switch (j) {
      case 0:
        return sin_taylor(t);
      case 1:
        return cos_taylor(t);
      case -1:
        return negative(cos_taylor(t));
      default:
        return negative(sin_taylor(t));
    }
  }

  GPU_dd u = d_dd_cos_table[abs_k-1];
  GPU_dd v = d_dd_sin_table[abs_k-1];
  GPU_dd sin_t, cos_t;
  sincos_taylor(t, sin_t, cos_t);
  if (j == 0) {
    if (k > 0) {
      r = u * sin_t + v * cos_t;
    } else {
      r = u * sin_t - v * cos_t;
    }
  } else if (j == 1) {
    if (k > 0) {
      r = u * cos_t - v * sin_t;
    } else {
      r = u * cos_t + v * sin_t;
    }
  } else if (j == -1) {
    if (k > 0) {
      r = v * sin_t - u * cos_t;
    } else if (k < 0) {
      //r = -u * cos_t - v * sin_t;
      r = negative(u * cos_t) - v * sin_t;
    }
  } else {
    if (k > 0) {
      //r = -u * sin_t - v * cos_t;
      r = negative(u * sin_t) - v * cos_t;
    } else {
      r = v * cos_t - u * sin_t;
    }
  }

  return r;
}

__device__
GPU_dd cos(const GPU_dd &a, const GPU_dd* d_dd_sin_table, const GPU_dd* d_dd_cos_table) {

  if (is_zero(a)) {
    return make_dd(1.0);
  }

  // approximately reduce modulo 2*pi
  GPU_dd z = nint(a / _dd_2pi);
  GPU_dd r = a - z * _dd_2pi;

  // approximately reduce modulo pi/2 and then modulo pi/16
  GPU_dd t;
  double q = floor(r.x / _dd_pi2.x + 0.5);
  t = r - _dd_pi2 * q;
  int j = (int)(q);
  q = floor(t.x / _dd_pi16.x + 0.5);
  t = t - _dd_pi16 * q;
  int k = (int)(q);
  int abs_k = abs(k);

  if (j < -2 || j > 2) {
    //dd_real::error("(dd_real::cos): Cannot reduce modulo pi/2.");
    //return dd_real::_nan;
    return make_dd(0.0);
  }

  if (abs_k > 4) {
    //dd_real::error("(dd_real::cos): Cannot reduce modulo pi/16.");
    //return dd_real::_nan;
    return make_dd(0.0);
  }

  if (k == 0) {
    switch (j) {
      case 0:
        return cos_taylor(t);
      case 1:
        return negative(sin_taylor(t));
      case -1:
        return sin_taylor(t);
      default:
        return negative(cos_taylor(t));
    }
  }

  GPU_dd sin_t, cos_t;
  sincos_taylor(t, sin_t, cos_t);
  GPU_dd u = d_dd_cos_table[abs_k-1];
  GPU_dd v = d_dd_sin_table[abs_k-1];

  if (j == 0) {
    if (k > 0) {
      r = u * cos_t - v * sin_t;
    } else {
      r = u * cos_t + v * sin_t;
    }
  } else if (j == 1) {
    if (k > 0) {
      r = negative(u * sin_t) - v * cos_t;
    } else {
      r = v * cos_t - u * sin_t;
    }
  } else if (j == -1) {
    if (k > 0) {
      r = u * sin_t + v * cos_t;
    } else {
      r = u * sin_t - v * cos_t;
    }
  } else {
    if (k > 0) {
      r = v * sin_t - u * cos_t;
    } else {
      r = negative(u * cos_t) - v * sin_t;
    }
  }

  return r;
}


#endif

