#include <iostream>
#include <math.h>

#define _fpu_getcw(cw) asm ("fnstcw %0" : "=m" (cw))
#define _fpu_setcw(cw) asm ("fldcw %0"  : "=m" (cw))

unsigned short int old_cw,new_cw;

// set the fpu to round down (towards -infty) 
#define _fpu_rndd {\
 _fpu_getcw(old_cw);\
 printf("old control word was %hx\n",old_cw);\
 new_cw=old_cw&0xF3FF;\
 new_cw=new_cw|0x0400;\
 _fpu_setcw(new_cw);}

#define _fpu_rndu {\
 _fpu_getcw(old_cw);\
 printf("old control word was %hx\n",old_cw);\
 new_cw=old_cw&0xF3FF;\
 new_cw=new_cw|0x0800;\
 _fpu_setcw(new_cw);}


#define _fpu_rndn {\
 _fpu_getcw(old_cw);\ 
 printf("old control word was %hx\n",old_cw);\
new_cw=old_cw&0xF3FF;\
 _fpu_setcw(new_cw);}

// restore the fpu control register from memory
#define _fpu_restore {_fpu_setcw(old_cw);}


using namespace std;

int main()
{
    double x,y;

    cout << "Enter x,y: ";
    cin >> x>> y;

  _fpu_rndn;
  printf("Default rounding.\n");
  printf("1/%20.18e is %20.18e\n",x,1.0/x);
  printf("sin(3.125 is %20.18e\n",sin(y));
  _fpu_rndd;
  printf("Rounding down.\n");
  printf("1/%20.18e is %20.18e\n",x,1.0/x);
  printf("sin(3.125 is %20.18e\n",sin(y));
  _fpu_rndu;
  printf("Rounding up.\n");
  printf("1/%20.18e is %20.18e\n",x,1.0/x);
  printf("sin(3.125 is %20.18e\n",sin(y));
  _fpu_rndn;
  return(0);
};
