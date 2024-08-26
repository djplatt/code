/* TFmode long double library function call
   interface to real.c  */


typedef		float TFtype	__attribute__ ((mode (TF)));

typedef struct {
  unsigned short e[10];
  } etype;
#define ETYPE etype

/* real.c references */
int ecmp ();
void eadd (), esub (), emul (), ediv();
void ltoe (), ultoe ();
void eifrac (), euifrac ();
void e24toe (), etoe24 ();
void e53toe (), etoe53 ();
void e113toe (), etoe113 ();

#define GETLONG e113toe
#define PUTLONG etoe113

TFtype
__addtf3 (TFtype a1, TFtype a2)
{
  ETYPE e1;
  ETYPE e2;
  ETYPE e3;
  TFtype a3;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  eadd (&e1, &e2, &e3);
  PUTLONG (&e3, &a3);
  return a3;
}

long double doadd( long double x, long double y )
{
long double z;

z = x + y;
return(z);
}


/* a1 - a2 */
TFtype
__subtf3 (TFtype a1, TFtype a2)
{
  ETYPE e1, e2, e3;
  TFtype a3;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  esub (&e2, &e1, &e3);
  PUTLONG (&e3, &a3);
  return a3;
}

long double
__multf3 (long double a1, long double a2)
{
  ETYPE e1, e2, e3;
  long double a3;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  emul (&e1, &e2, &e3);
  PUTLONG (&e3, &a3);
  return a3;
}

/* a1 divided by a2 */
long double
__divtf3 (long double a1, long double a2)
{
  ETYPE e1, e2, e3;
  long double a3;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  ediv (&e2, &e1, &e3);
  PUTLONG (&e3, &a3);
  return a3;
}


long double
__negtf2 (long double a1)
{
  ETYPE e1;
  long double a2;

  GETLONG (&a1, &e1);
  eneg (&e1);
  PUTLONG (&e1, &a2);
  return a2;
}


/* compare two floats
 * return -1 if a1 < a2
 *        +1 if a1 > a2
 *         0 if a1 == a2
 */
long
__cmptf2 (long double a1, long double a2)
{
  ETYPE e1, e2;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  return ((long )ecmp(&e1, &e2));
}


int
__eqtf2 (long double a1, long double a2)
{
  ETYPE e1, e2;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  return (ecmp(&e1, &e2) == 0);
}


int
__netf2 (long double a1, long double a2)
{
  ETYPE e1, e2;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  return (ecmp(&e1, &e2) != 0);
}


int
__gttf2 (long double a1, long double a2)
{
  ETYPE e1, e2;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  return (ecmp(&e1, &e2) == 1);
}


int
__getf2 (long double a1, long double a2)
{
  ETYPE e1, e2;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  return (ecmp(&e1, &e2) >= 0);
}


int
__lttf2 (long double a1, long double a2)
{
  ETYPE e1, e2;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  return (ecmp(&e1, &e2) == -1);
}


int
__letf2 (long double a1, long double a2)
{
  ETYPE e1, e2;
  int k;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  k = ecmp (&e1, &e2);
  return ((k == 0) || (k == -1));
}


long double
__floatsitf (long a1)
{
  ETYPE e1;
  long double a2;

  ltoe (&a1, &e1);
  PUTLONG (&e1, &a2);
  return a2;
}


long double
__floatunssitf (unsigned long a1)
{
  ETYPE e1;
  long double a2;

  ultoe (&a1, &e1);
  PUTLONG (&e1, &a2);
  return a2;
}


unsigned long
__fixunstfsi (long double a1)
{
  ETYPE e1, e2;
  unsigned long ul;

  GETLONG (&a1, &e1);
  euifrac (&e1, &ul, &e2);
  return ul;
}

/*
eifrac( e, &l, frac )   e to long integer and e type fraction
euifrac( e, &l, frac )  e to unsigned long integer and e type fraction
*/

long
__fixtfsi (long double a1)
{
  ETYPE e1, e2;
  unsigned long ul;

  GETLONG (&a1, &e1);
  eifrac (&e1, &ul, &e2);
  return ul;
}

/* convert long double to float */
float
__trunctfsf2 (long double a1)
{
  ETYPE e1;
  float f1;

  GETLONG (&a1, &e1);
  etoe24 (&e1, &f1);
  return f1;
}

/* convert long double to float */
double
__trunctfdf2 (long double a1)
{
  ETYPE e1;
  double d1;

  GETLONG (&a1, &e1);
  etoe53 (&e1, &d1);
  return d1;
}


long double
__extendsftf2 (float a1)
{
  ETYPE e1;
  long double a2;

  e24toe (&a1, &e1);
  PUTLONG (&e1, &a2);
  return a2;
}


long double
__extenddftf2 (double a1)
{
  ETYPE e1;
  long double a2;

  e53toe (&a1, &e1);
  PUTLONG (&e1, &a2);
  return a2;
}
