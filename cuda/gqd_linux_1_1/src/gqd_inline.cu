#ifndef _GQD_INLINE_CU_
#define _GQD_INLINE_CU_



/** normalization functions */
__host__ __device__
void quick_renorm(double &c0, double &c1, 
				  double &c2, double &c3, double &c4) 
{
	double t0, t1, t2, t3;
	double s;
	s  = quick_two_sum(c3, c4, t3);
	s  = quick_two_sum(c2, s , t2);
	s  = quick_two_sum(c1, s , t1);
	c0 = quick_two_sum(c0, s , t0);

	s  = quick_two_sum(t2, t3, t2);
	s  = quick_two_sum(t1, s , t1);
	c1 = quick_two_sum(t0, s , t0);

	s  = quick_two_sum(t1, t2, t1);
	c2 = quick_two_sum(t0, s , t0);

	c3 = t0 + t1;
}

__host__ __device__
void renorm(double &c0, double &c1, 
			double &c2, double &c3) 
{
	double s0, s1, s2 = 0.0, s3 = 0.0;

	//if (isinf(c0)) return;

	s0 = quick_two_sum(c2, c3, c3);
	s0 = quick_two_sum(c1, s0, c2);
	c0 = quick_two_sum(c0, s0, c1);

	s0 = c0;
	s1 = c1;
	if (s1 != 0.0) {
		s1 = quick_two_sum(s1, c2, s2);
		if (s2 != 0.0)
			s2 = quick_two_sum(s2, c3, s3);
		else
			s1 = quick_two_sum(s1, c3, s2);
	} else {
		s0 = quick_two_sum(s0, c2, s1);
		if (s1 != 0.0)
			s1 = quick_two_sum(s1, c3, s2);
		else
			s0 = quick_two_sum(s0, c3, s1);
	}

	c0 = s0;
	c1 = s1;
	c2 = s2;
	c3 = s3;
}

__host__ __device__
void renorm(double &c0, double &c1, 
			double &c2, double &c3, double &c4) 
{
	double s0, s1, s2 = 0.0, s3 = 0.0;

	//if (isinf(c0)) return;

	s0 = quick_two_sum(c3, c4, c4);
	s0 = quick_two_sum(c2, s0, c3);
	s0 = quick_two_sum(c1, s0, c2);
	c0 = quick_two_sum(c0, s0, c1);

	s0 = c0;
	s1 = c1;

	s0 = quick_two_sum(c0, c1, s1);
	if (s1 != 0.0) 
	{
		s1 = quick_two_sum(s1, c2, s2);
		if (s2 != 0.0) {
			s2 = quick_two_sum(s2, c3, s3);
			if (s3 != 0.0)
				s3 += c4;
			else
				s2 += c4;
		} else {
			s1 = quick_two_sum(s1, c3, s2);
			if (s2 != 0.0)
				s2 = quick_two_sum(s2, c4, s3);
			else
				s1 = quick_two_sum(s1, c4, s2);
		}
	} else {
		s0 = quick_two_sum(s0, c2, s1);
		if (s1 != 0.0) {
			s1 = quick_two_sum(s1, c3, s2);
			if (s2 != 0.0)
				s2 = quick_two_sum(s2, c4, s3);
			else
				s1 = quick_two_sum(s1, c4, s2);
		} else {
			s0 = quick_two_sum(s0, c3, s1);
			if (s1 != 0.0)
				s1 = quick_two_sum(s1, c4, s2);
			else
				s0 = quick_two_sum(s0, c4, s1);
		}
	}

	c0 = s0;
	c1 = s1;
	c2 = s2;
	c3 = s3;
}

__host__ __device__
void renorm( GPU_qd &x ) {
	renorm(x.d1.x, x.d1.y, x.d2.x, x.d2.y);
}

__host__ __device__
void renorm( GPU_qd &x, double &e) {
	renorm(x.d1.x, x.d1.y, x.d2.x, x.d2.y, e);
}

/** additions */
__host__ __device__
void three_sum(double &a, double &b, double &c)
{
        double t1, t2, t3;
        t1 = two_sum(a, b, t2);
        a  = two_sum(c, t1, t3);
        b  = two_sum(t2, t3, c);
}

__host__ __device__
void three_sum2(double &a, double &b, double &c) {
        double t1, t2, t3;
        t1 = two_sum(a, b, t2);
        a  = two_sum(c, t1, t3);
        b = (t2 + t3);
}

///qd = qd + double
__host__ __device__
GPU_qd operator+(const GPU_qd &a, double b) {
        double c0, c1, c2, c3;
        double e;

        c0 = two_sum(a.d1.x, b, e);
        c1 = two_sum(a.d1.y, e, e);
        c2 = two_sum(a.d2.x, e, e);
        c3 = two_sum(a.d2.y, e, e);

        renorm(c0, c1, c2, c3, e);

        return make_qd(c0, c1, c2, c3);
}

///qd = double + qd
__host__ __device__
GPU_qd operator+( double a, const GPU_qd &b )
{
        return ( b + a );
}

///qd = qd + qd
__host__ __device__
GPU_qd sloppy_add(const GPU_qd &a, const GPU_qd &b)
{
        double s0, s1, s2, s3;
        double t0, t1, t2, t3;

        double v0, v1, v2, v3;
        double u0, u1, u2, u3;
        double w0, w1, w2, w3;

        s0 = a.d1.x + b.d1.x;
        s1 = a.d1.y + b.d1.y;
        s2 = a.d2.x + b.d2.x;
        s3 = a.d2.y + b.d2.y;

        v0 = s0 - a.d1.x;
        v1 = s1 - a.d1.y;
        v2 = s2 - a.d2.x;
        v3 = s3 - a.d2.y;

        u0 = s0 - v0;
        u1 = s1 - v1;
        u2 = s2 - v2;
        u3 = s3 - v3;

        w0 = a.d1.x - u0;
        w1 = a.d1.y - u1;
        w2 = a.d2.x - u2;
        w3 = a.d2.y - u3;

        u0 = b.d1.x - v0;
        u1 = b.d1.y - v1;
        u2 = b.d2.x - v2;
        u3 = b.d2.y - v3;

        t0 = w0 + u0;
        t1 = w1 + u1;
        t2 = w2 + u2;
        t3 = w3 + u3;

        s1 = two_sum(s1, t0, t0);
        three_sum(s2, t0, t1);
        three_sum2(s3, t0, t2);
        t0 = t0 + t1 + t3;
  
	renorm(s0, s1, s2, s3, t0);

        return make_qd(s0, s1, s2, s3);
}

__host__ __device__
GPU_qd operator+(const GPU_qd &a, const GPU_qd &b)
{
        return sloppy_add(a, b);
}


/** subtractions */
__host__ __device__
GPU_qd negative( const GPU_qd &a )
{
        return make_qd( -a.d1.x, -a.d1.y, -a.d2.x, -a.d2.y );
}

__host__ __device__
GPU_qd operator-(const GPU_qd &a, double b)
{
        return (a + (-b));
}

__host__ __device__
GPU_qd operator-(double a, const GPU_qd &b)
{
        return (a + negative(b));
}

__host__ __device__
GPU_qd operator-(const GPU_qd &a, const GPU_qd &b)
{
        return (a + negative(b));
}

/** multiplications */
__host__ __device__
GPU_qd mul_pwr2(const GPU_qd &a, double b) {
        return make_qd(a.d1.x * b, a.d1.y * b, a.d2.x * b, a.d2.y * b);
}


//quad_double * double
 __device__
GPU_qd operator*(const GPU_qd &a, double b)
{
        double p0, p1, p2, p3;
        double q0, q1, q2;
        double s0, s1, s2, s3, s4;

        p0 = two_prod(a.d1.x, b, q0);
        p1 = two_prod(a.d1.y, b, q1);
        p2 = two_prod(a.d2.x, b, q2);
        p3 = a.d2.y * b;

        s0 = p0;

        s1 = two_sum(q0, p1, s2);

        three_sum(s2, q1, p2);

        three_sum2(q1, q2, p3);
        s3 = q1;

        s4 = q2 + p2;

        renorm(s0, s1, s2, s3, s4);
        return make_qd(s0, s1, s2, s3);
}
//quad_double = double*quad_double
__device__
GPU_qd operator*( double a, const GPU_qd &b )
{
        return b*a;
}

__device__
GPU_qd sloppy_mul(const GPU_qd &a, const GPU_qd &b)
{
        double p0, p1, p2, p3, p4, p5;
        double q0, q1, q2, q3, q4, q5;
        double t0, t1;
        double s0, s1, s2;

        p0 = two_prod(a.d1.x, b.d1.x, q0);

        p1 = two_prod(a.d1.x, b.d1.y, q1);
        p2 = two_prod(a.d1.y, b.d1.x, q2);

        p3 = two_prod(a.d1.x, b.d2.x, q3);
        p4 = two_prod(a.d1.y, b.d1.y, q4);
        p5 = two_prod(a.d2.x, b.d1.x, q5);


        /* Start Accumulation */
        three_sum(p1, p2, q0);

	//return make_qd(p1, p2, q0, 0.0);

        /* Six-Three Sum  of p2, q1, q2, p3, p4, p5. */
        three_sum(p2, q1, q2);
        three_sum(p3, p4, p5);
        /* compute (s0, s1, s2) = (p2, q1, q2) + (p3, p4, p5). */
        s0 = two_sum(p2, p3, t0);
        s1 = two_sum(q1, p4, t1);
	s2 = q2 + p5;
        s1 = two_sum(s1, t0, t0);
        s2 += (t0 + t1);

	//return make_qd(s0, s1, t0, t1);

        /* O(eps^3) order terms */
        //!!!s1 = s1 + (a.d1.x*b.d2.y + a.d1.y*b.d2.x + a.d2.x*b.d1.y + a.d2.y*b.d1.x + q0 + q3 + q4 + q5);
	
	s1 = s1 + (__dmul_rn(a.d1.x,b.d2.y) + __dmul_rn(a.d1.y,b.d2.x) + 
			__dmul_rn(a.d2.x,b.d1.y) + __dmul_rn(a.d2.y,b.d1.x) + q0 + q3 + q4 + q5);
	renorm(p0, p1, s0, s1, s2);

        return make_qd(p0, p1, s0, s1);
	
}

 __device__
GPU_qd operator*(const GPU_qd &a, const GPU_qd &b) {
        return sloppy_mul(a, b);
}

 __device__
GPU_qd sqr(const GPU_qd &a) 
{
	double p0, p1, p2, p3, p4, p5;
	double q0, q1, q2, q3;
	double s0, s1;
	double t0, t1;

	p0 = two_sqr(a.d1.x, q0);
	p1 = two_prod(2.0 * a.d1.x, a.d1.y, q1);
	p2 = two_prod(2.0 * a.d1.x, a.d2.x, q2);
	p3 = two_sqr(a.d1.y, q3);

	p1 = two_sum(q0, p1, q0);

	q0 = two_sum(q0, q1, q1);
	p2 = two_sum(p2, p3, p3);

	s0 = two_sum(q0, p2, t0);
	s1 = two_sum(q1, p3, t1);

	s1 = two_sum(s1, t0, t0);
	t0 += t1;

	s1 = quick_two_sum(s1, t0, t0);
	p2 = quick_two_sum(s0, s1, t1);
	p3 = quick_two_sum(t1, t0, q0);

	p4 = 2.0 * a.d1.x * a.d2.y;
	p5 = 2.0 * a.d1.y * a.d2.x;

	p4 = two_sum(p4, p5, p5);
	q2 = two_sum(q2, q3, q3);

	t0 = two_sum(p4, q2, t1);
	t1 = t1 + p5 + q3;

	p3 = two_sum(p3, t0, p4);
	p4 = p4 + q0 + t1;

	renorm(p0, p1, p2, p3, p4);
	return make_qd(p0, p1, p2, p3);
}

/** divisions */
__device__
GPU_qd sloppy_div(const GPU_qd &a, const GPU_qd &b) 
{
	double q0, q1, q2, q3;

	GPU_qd r;

	q0 = a.d1.x / b.d1.x;
	r = a - (b * q0);

	q1 = r.d1.x / b.d1.x;
	r = r - (b * q1);

	q2 = r.d1.x / b.d1.x;
	r = r - (b * q2);

	q3 = r.d1.x / b.d1.x;

	renorm(q0, q1, q2, q3);

	return make_qd(q0, q1, q2, q3);
}

__device__
GPU_qd operator/(const GPU_qd &a, const GPU_qd &b) 
{
	return sloppy_div(a, b);
}

/* double / quad-double */
__device__
GPU_qd operator/(double a, const GPU_qd &b) 
{
	return make_qd(a) / b;
}

/* quad-double / double */
__device__
GPU_qd operator/( const GPU_qd &a, double b )
{
	return a/make_qd(b);
}

/********** Miscellaneous **********/
__host__ __device__
GPU_qd abs(const GPU_qd &a) 
{
	return (a.d1.x < 0.0) ? (negative(a)) : (a);
}

/********************** Simple Conversion ********************/
__host__ __device__
double to_double(const GPU_qd &a) 
{
	return a.d1.x;
}

__host__ __device__
GPU_qd ldexp(const GPU_qd &a, int n) 
{
	return make_qd(ldexp(a.d1.x, n), ldexp(a.d1.y, n), 
		ldexp(a.d2.x, n), ldexp(a.d2.y, n));
}

__device__
GPU_qd inv(const GPU_qd &qd) 
{
	return 1.0 / qd;
}


/********** Greater-Than Comparison ***********/

__host__ __device__
bool operator>=(const GPU_qd &a, const GPU_qd &b) 
{
	return (a.d1.x > b.d1.x || 
		(a.d1.x == b.d1.x && (a.d1.y > b.d1.y ||
		(a.d1.y == b.d1.y && (a.d2.x > b.d2.x ||
		(a.d2.x == b.d2.x && a.d2.y >= b.d2.y))))));
}


/********** Less-Than Comparison ***********/
__host__ __device__
bool operator<(const GPU_qd &a, double b) {
	return (a.d1.x < b || (a.d1.x == b && a.d1.y < 0.0));
}

__host__ __device__
bool operator<(const GPU_qd &a, const GPU_qd &b) {
	return (a.d1.x < b.d1.x ||
		(a.d1.x == b.d1.x && (a.d1.y < b.d1.y ||
		(a.d1.y == b.d1.y && (a.d2.x < b.d2.x ||
		(a.d2.x == b.d2.x && a.d2.y < b.d2.y))))));
}

__host__ __device__
bool operator<=(const GPU_qd &a, const GPU_qd &b) {
  return (a.d1.x < b.d1.x || 
          (a.d1.x == b.d1.x && (a.d1.y < b.d1.y ||
                            (a.d1.y == b.d1.y && (a.d2.x < b.d2.x ||
                                              (a.d2.x == b.d2.x && a.d2.y <= b.d2.y))))));
}

__host__ __device__
bool operator==(const GPU_qd &a, const GPU_qd &b) {
  return (a.d1.x == b.d1.x && a.d1.y == b.d1.y && 
          a.d2.x == b.d2.x && a.d2.y == b.d2.y);
}



/********** Greater-Than Comparison ***********/
__host__ __device__
bool operator>(const GPU_qd &a, double b) {
	return (a.d1.x > b || (a.d1.x == b && a.d1.y > 0.0));
}

__host__ __device__
bool operator<(double a, const GPU_qd &b) {
	return (b > a);
}

__host__ __device__
bool operator>(double a, const GPU_qd &b) {
	return (b < a);
}

__host__ __device__ 
bool operator>(const GPU_qd &a, const GPU_qd &b) {
	return (a.d1.x > b.d1.x ||
		(a.d1.x == b.d1.x && (a.d1.y > b.d1.y ||
		(a.d1.y == b.d1.y && (a.d2.x > b.d2.x ||
		(a.d2.x == b.d2.x && a.d2.y > b.d2.y))))));
}


__host__ __device__
bool is_zero( const GPU_qd &x ) 
{
	return (x.d1.x == 0.0);
}

__host__ __device__
bool is_one( const GPU_qd &x ) 
{
	return (x.d1.x == 1.0 && x.d1.y == 0.0 && x.d2.x == 0.0 && x.d2.y == 0.0);
}

__host__ __device__
bool is_positive( const GPU_qd &x )
{
	return (x.d1.x > 0.0);
}

__host__ __device__
bool is_negative( const GPU_qd &x ) 
{
	return (x.d1.x < 0.0);
}

__device__
GPU_qd nint(const GPU_qd &a) {
  double x0, x1, x2, x3;

  x0 = nint(a.d1.x);
  x1 = x2 = x3 = 0.0;

  if (x0 == a.d1.x) {
    /* First double is already an integer. */
    x1 = nint(a.d1.y);

    if (x1 == a.d1.y) {
      /* Second double is already an integer. */
      x2 = nint(a.d2.x);
      
      if (x2 == a.d2.x) {
        /* Third double is already an integer. */
        x3 = nint(a.d2.y);
      } else {
        if (fabs(x2 - a.d2.x) == 0.5 && a.d2.y < 0.0) {
          x2 -= 1.0;
        }
      }

    } else {
      if (fabs(x1 - a.d1.y) == 0.5 && a.d2.x < 0.0) {
          x1 -= 1.0;
      }
    }

  } else {
    /* First double is not an integer. */
      if (fabs(x0 - a.d1.x) == 0.5 && a.d1.y < 0.0) {
          x0 -= 1.0;
      }
  }
  
  renorm(x0, x1, x2, x3);
  return make_qd(x0, x1, x2, x3);
}

#endif


