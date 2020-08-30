INCLUDE=/data/local/include
LIBRARY=/data/local/lib
all: dirichlet dirichlet_double

dirichlet:	dirichlet1.0.cpp 
	g++ dirichlet1.0.cpp -o dirichlet -finline-functions -fomit-frame-pointer -I ${INCLUDE} -L ${LIBRARY} -larb -lflint -lgmp

dirichlet_double:	dirichlet_double.cpp 
	g++ dirichlet_double.cpp -o dirichlet_double -finline-functions -fomit-frame-pointer -I ${INCLUDE} -L ${LIBRARY} -larb -lflint -lgmp

l-func-old:	main8.0.cpp structs9.0.h error9.0.h acb_fft.h
	g++ main8.0.cpp -o l-func-old -finline-functions -fomit-frame-pointer -I ${INCLUDE} -L ${LIBRARY} -larb -lflint -lgmp


l-func-hi:	main10.0_hi.cpp main10.0.cpp structs9.0.h error9.0.h acb_fft.h
	g++ main10.0_hi.cpp -o l-func-hi -finline-functions -fomit-frame-pointer -I ${INCLUDE} -L ${LIBRARY} -larb -lflint -lgmp
