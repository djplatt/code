all:	g_sum g_sum1 A_T first_N hassani hassani_arb

g_sum:	g_sum.cpp
	g++ -finline-functions -fomit-frame-pointer -frounding-math -mfpmath=387 -o g_sum g_sum.cpp -lcrlibm

g_sum1:	g_sum1.cpp
	g++ -finline-functions -fomit-frame-pointer -frounding-math -mfpmath=387 -o g_sum1 g_sum1.cpp -lcrlibm

A_T:	A_T.cpp
	g++ -finline-functions -fomit-frame-pointer -frounding-math -mfpmath=387 -o A_T A_T.cpp -lcrlibm


first_N:	first_N.cpp
	g++ -finline-functions -fomit-frame-pointer -frounding-math -mfpmath=387 -o first_N first_N.cpp -lcrlibm

hassani:	hassani.cpp
	g++ -finline-functions -fomit-frame-pointer -frounding-math -mfpmath=387 -o hassani hassani.cpp -lcrlibm

hassani_arb:	hassani_arb.c
	gcc -O2 hassani_arb.c -o hassani_arb -larb