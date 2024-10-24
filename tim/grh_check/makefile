all:	dir_even_deriv dir_odd_deriv

dir_even_deriv:	dir_even_deriv.c
	gcc -O2 -o dir_even_deriv dir_even_deriv.c -I ~/glfunc_stripdown/include/ -L ~/glfunc_stripdown/build/ -llfun -lflint

dir_odd_deriv:	dir_odd_deriv.c
	gcc -O2 -o dir_odd_deriv dir_odd_deriv.c -I ~/glfunc_stripdown/include/ -L ~/glfunc_stripdown/build/ -llfun -lflint
