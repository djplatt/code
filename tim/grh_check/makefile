all:	dir_even_deriv dir_odd_deriv

dir_even_deriv:	dir_even_deriv.c
	gcc -O2 -o dir_even_deriv dir_even_deriv.c -I ~/glfunc_stripdown/include/ -I ~/local/include ~/glfunc_stripdown/build/liblfun.so ~/local/lib/libflint.so

dir_odd_deriv:	dir_odd_deriv.c
	gcc -O2 -o dir_odd_deriv dir_odd_deriv.c -I ~/glfunc_stripdown/include/ -I ~/local/include ~/glfunc_stripdown/build/liblfun.so ~/local/lib/libflint.so
