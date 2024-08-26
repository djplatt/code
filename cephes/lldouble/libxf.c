/* XFmode long double library function call
   interface to real.c  */



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

long double
__addxf3 (long double a1, long double a2)
{
  ETYPE e1, e2, e3;
  long double a3;

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
long double
__subxf3 (long double a1, long double a2)
{
  ETYPE e1, e2, e3;
  long double a3;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  esub (&e2, &e1, &e3);
  PUTLONG (&e3, &a3);
  return a3;
}

long double
__mulxf3 (long double a1, long double a2)
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
__divxf3 (long double a1, long double a2)
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
__negxf2 (long double a1)
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
__cmpxf2 (long double a1, long double a2)
{
  ETYPE e1, e2;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  return ((long )ecmp(&e1, &e2));
}


int
__eqxf2 (long double a1, long double a2)
{
  ETYPE e1, e2;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  return (ecmp(&e1, &e2) == 0);
}


int
__nexf2 (long double a1, long double a2)
{
  ETYPE e1, e2;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  return (ecmp(&e1, &e2) != 0);
}


int
__gtxf2 (long double a1, long double a2)
{
  ETYPE e1, e2;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  return (ecmp(&e1, &e2) == 1);
}


int
__gexf2 (long double a1, long double a2)
{
  ETYPE e1, e2;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  return (ecmp(&e1, &e2) >= 0);
}


int
__ltxf2 (long double a1, long double a2)
{
  ETYPE e1, e2;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  return (ecmp(&e1, &e2) == -1);
}


int
__lexf2 (long double a1, long double a2)
{
  ETYPE e1, e2;
  int k;

  GETLONG (&a1, &e1);
  GETLONG (&a2, &e2);
  k = ecmp (&e1, &e2);
  return ((k == 0) || (k == -1));
}


long double
__floatsixf (long a1)
{
  ETYPE e1;
  long double a2;

  ltoe (&a1, &e1);
  PUTLONG (&e1, &a2);
  return a2;
}


long double
__floatunssixf (unsigned long a1)
{
  ETYPE e1;
  long double a2;

  ultoe (&a1, &e1);
  PUTLONG (&e1, &a2);
  return a2;
}


unsigned long
__fixunsxfsi (long double a1)
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
__fixxfsi (long double a1)
{
  ETYPE e1, e2;
  unsigned long ul;

  GETLONG (&a1, &e1);
  eifrac (&e1, &ul, &e2);
  return ul;
}

/* convert long double to float */
float
__truncxfsf2 (long double a1)
{
  ETYPE e1;
  float f1;

  GETLONG (&a1, &e1);
  etoe24 (&e1, &f1);
  return f1;
}

/* convert long double to float */
double
__truncxfdf2 (long double a1)
{
  ETYPE e1;
  double d1;

  GETLONG (&a1, &e1);
  etoe53 (&e1, &d1);
  return d1;
}


long double
__extendsfxf2 (float a1)
{
  ETYPE e1;
  long double a2;

  e24toe (&a1, &e1);
  PUTLONG (&e1, &a2);
  return a2;
}


long double
__extenddfxf2 (double a1)
{
  ETYPE e1;
  long double a2;

  e53toe (&a1, &e1);
  PUTLONG (&e1, &a2);
  return a2;
}
