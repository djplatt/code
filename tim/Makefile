CPPOPTS=-O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387
LPPOPTS= -lcrlibm -lprimesieve -lmpfi
CEXES=
CPPEXES=cheb fitq lam1.0 tim3.0 zetalim3.0
EXES=$(CEXES) $(CPPEXES)
CPPSRCS=cheb.cpp fitq.cpp lam1.0.cpp tim3.0.cpp zetalim3.0.cpp
CPPOBJS=$(subst .cpp,.o,$(CPPSRCS))
CPPINCLUDES=../includes/int_double12.0.h

all: $(EXES)

$(CPPEXES): % : %.o
	g++ $< -o $@ $(LPPOPTS) 

$(CPPOBJS): %.o : %.cpp $(CPPINCLUDES)
	g++ $(CPPOPTS) -c $< -o $@

clean:
	\rm -f $(CPPOBJS) $(COBJS)

spotless:
	\rm -f $(EXES)
	make clean
