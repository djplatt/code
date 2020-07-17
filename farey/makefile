all: farey_int farey_threads farey_threads_int

farey_int: farey_int.cpp ../includes/int_double13.0.h
	g++ farey_int.cpp -o farey_int -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ~/crlibm-1.0beta4/ -L ~/crlibm-1.0beta4/ -lcrlibm -lm

farey_threads_int: farey_threads_int.cpp ../includes/int_double13.0.h
	g++ farey_threads_int.cpp -o farey_threads_int -O1 -pthread -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ~/crlibm-1.0beta4/ -L ~/crlibm-1.0beta4/ -lcrlibm -lm

farey_threads: farey_threads.cpp slint.h
	g++ -O3 -o farey_threads -pthread farey_threads.cpp

clean:
	rm -f farey_int
	rm -f farey_threads
	rm -f farey_threads_int

