all:	S1 S1_v1 S168_pi

S1:	S1.c
	gcc -O2 -I ${HOME}/code/quad -o S1 S1.c -larb -L ${HOME}/code/quad -lquad

S168_pi:	S168_pi.c
	gcc -O2 -I ${HOME}/code/quad -o S168_pi S168_pi.c -larb -L ${HOME}/code/quad -lquad

S1_v1:	S1_v1.c
	gcc -O2 -I ${HOME}/code/quad -o S1_v1 S1_v1.c -larb -L ${HOME}/code/quad -lquad

clean:
	rm -f S1 S1_v1
