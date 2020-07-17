all: zero_stats debug start_end sum more_sums buethe_check zero_gaps artur artur_sums artur_ints artur_arb_sum yet_more_sums

artur_sums: artur_sums.cpp
	g++ artur_sums.cpp -o artur_sums -O2 -finline-functions -fomit-frame-pointer -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.15.1/ -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.15.1 -larb

artur_arb_sum: artur_arb_sum.cpp
	g++ artur_arb_sum.cpp -o artur_arb_sum -O2 -finline-functions -fomit-frame-pointer -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.15.1/ -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.15.1 -larb

artur_ints: artur_ints.cpp
	g++ artur_ints.cpp -o artur_ints -O2 -finline-functions -fomit-frame-pointer -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.15.1/ -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.15.1 -larb


artur: artur1.0.cpp
	g++ artur1.0.cpp -o artur -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ${HOME}/crlibm-1.0beta4/ -L ~/crlibm-1.0beta4/ -lcrlibm -lm

zero_gaps: zero_gaps.cpp ../includes/int_double12.0.h
	g++ zero_gaps.cpp -o zero_gaps -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ~/crlibm-1.0beta4/ -L ~/crlibm-1.0beta4/ -lcrlibm -lm

zero_stats: zero_stats1.0.cpp ../includes/int_double12.0.h ../includes/hurwitz1.0.h
	g++ zero_stats1.0.cpp -o zero_stats -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ~/crlibm-1.0beta4/ -L ~/crlibm-1.0beta4/ -lcrlibm -lm

buethe_check: buethe_check.cpp ../includes/int_double12.0.h ../includes/hurwitz1.0.h
	g++ buethe_check.cpp -o buethe_check -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ~/crlibm-1.0beta4/ -L ~/crlibm-1.0beta4/ -lcrlibm -lm

debug: debug.cpp ../includes/int_double12.0.h ../includes/hurwitz1.0.h
	g++ debug.cpp -o debug -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ~/crlibm-1.0beta4/ -L ~/crlibm-1.0beta4/ -lcrlibm -lm


start_end: start_end.cpp ../includes/int_double12.0.h ../includes/hurwitz1.0.h
	g++ start_end.cpp -o start_end -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ~/crlibm-1.0beta4/ -L ~/crlibm-1.0beta4/ -lcrlibm -lm

sum: sum.cpp ../includes/int_double12.0.h
	g++ sum.cpp -o sum -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ~/crlibm-1.0beta4/ -L ~/crlibm-1.0beta4/ -lcrlibm -lm

more_sums: more_sums.cpp ../includes/int_double12.0.h
	g++ more_sums.cpp -o more_sums -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ~/crlibm-1.0beta4/ -L ~/crlibm-1.0beta4/ -lcrlibm -lm

yet_more_sums: yet_more_sums.cpp ../includes/int_double14.2.h
	g++ yet_more_sums.cpp -o yet_more_sums -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ~/crlibm-1.0beta4/ -L ~/crlibm-1.0beta4/ -lcrlibm -lm


clean:
	rm -f zero_stats  debug start_end sum more_sums artur_ints artur_arb_sum yet_more_sums