EXES=xi
SRCS=xi.c
MPFI_LIB=$(HOME)/mpfi-1.5.1/src/.libs
MPFI_INCLUDE=$(HOME)/mpfi-1.5.1/src
#OBJS=$(subst .c,.o,$(SRCS))
#INCLUDES=../includes/mpfi_c.h ../includes/mpfi_fft.h ../includes/co_prime.h win_zeta.h

#CRDIR=$(HOME)/crlibm-1.0beta4

#PARI_INCLUDE=/usr/local/apps/pari-2.4.3.alpha/include/pari
#PARI_LIB=/usr/local/apps/pari-2.4.3.alpha/lib

COPTS=-O1 -frounding-math -fomit-frame-pointer -finline-functions -I $(MPFI_INCLUDE)

LOPTS=-L $(MPFI_LIB) -lmpfi -lm

all: $(EXES)

# This is a 'static pattern rule'
$(EXES): % : %.c 
	gcc $(COPTS) $(LOPTS) $< -o $@


.PHONY: clean spotless

clean:
	\rm -f $(EXES)

