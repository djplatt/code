	.file	"debug.cpp"
	.text
.globl _Z5width10int_double
	.type	_Z5width10int_double, @function
_Z5width10int_double:
.LFB1056:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	movsd	%xmm0, -32(%rsp)
	movsd	%xmm1, -48(%rsp)
	fldl	-48(%rsp)
	fchs
	fsubl	-32(%rsp)
	fstpl	-32(%rsp)
	movsd	-32(%rsp), %xmm0
	ret
	.cfi_endproc
.LFE1056:
	.size	_Z5width10int_double, .-_Z5width10int_double
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC1:
	.string	"[ %20.18e , %20.18e ]"
	.text
.globl _Z16print_int_doubleRK10int_double
	.type	_Z16print_int_doubleRK10int_double, @function
_Z16print_int_doubleRK10int_double:
.LFB1055:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	subq	$24, %rsp
	.cfi_def_cfa_offset 32
	fldl	8(%rdi)
	fchs
	fstpl	8(%rsp)
	movsd	8(%rsp), %xmm1
	movsd	(%rdi), %xmm0
	movl	$.LC1, %edi
	movl	$2, %eax
	call	printf
	addq	$24, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE1055:
	.size	_Z16print_int_doubleRK10int_double, .-_Z16print_int_doubleRK10int_double
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC2:
	.string	"Error constructing int_double, right %20.18e < left %20.18e . Exiting.\n"
	.section	.text._ZN10int_doubleC2Edd,"axG",@progbits,_ZN10int_doubleC5Edd,comdat
	.align 2
	.weak	_ZN10int_doubleC2Edd
	.type	_ZN10int_doubleC2Edd, @function
_ZN10int_doubleC2Edd:
.LFB1053:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	subq	$24, %rsp
	.cfi_def_cfa_offset 32
	movsd	%xmm1, 8(%rsp)
	fldl	8(%rsp)
	movsd	%xmm0, 8(%rsp)
	fldl	8(%rsp)
	fucomip	%st(1), %st
	jbe	.L10
	movapd	%xmm0, %xmm1
	fstpl	8(%rsp)
	movsd	8(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L10:
	movsd	%xmm0, (%rdi)
	fchs
	fstpl	8(%rdi)
	addq	$24, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE1053:
	.size	_ZN10int_doubleC2Edd, .-_ZN10int_doubleC2Edd
	.weak	_ZN10int_doubleC1Edd
	.set	_ZN10int_doubleC1Edd,_ZN10int_doubleC2Edd
	.text
	.type	_GLOBAL__I__Z16print_int_doubleRK10int_double, @function
_GLOBAL__I__Z16print_int_doubleRK10int_double:
.LFB1191:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	subq	$88, %rsp
	.cfi_def_cfa_offset 96
	movl	$_ZStL8__ioinit, %edi
	call	_ZNSt8ios_base4InitC1Ev
	movl	$__dso_handle, %edx
	movl	$_ZStL8__ioinit, %esi
	movl	$_ZNSt8ios_base4InitD1Ev, %edi
	call	__cxa_atexit
	flds	.LC3(%rip)
	fstpl	8(%rsp)
	movsd	8(%rsp), %xmm1
	movapd	%xmm1, %xmm0
	movl	$_ZL6d_half, %edi
	call	_ZN10int_doubleC1Edd
	movq	_d_pi2(%rip), %rax
	movsd	(%rax), %xmm1
	movq	_d_pi(%rip), %rax
	movsd	(%rax), %xmm0
	movl	$d_pi, %edi
	call	_ZN10int_doubleC1Edd
	movq	_d_gamma2(%rip), %rax
	movsd	(%rax), %xmm1
	movq	_d_gamma(%rip), %rax
	movsd	(%rax), %xmm0
	movl	$d_gamma, %edi
	call	_ZN10int_doubleC1Edd
	movq	_d_pi2_2(%rip), %rax
	movsd	(%rax), %xmm1
	movq	_d_pi_2(%rip), %rax
	movsd	(%rax), %xmm0
	movl	$d_pi_2, %edi
	call	_ZN10int_doubleC1Edd
	movq	_d_2pi2_2(%rip), %rax
	movsd	(%rax), %xmm1
	movq	_d_2pi_2(%rip), %rax
	movsd	(%rax), %xmm0
	movl	$d_two_pi, %edi
	call	_ZN10int_doubleC1Edd
	movq	__nze(%rip), %rax
	movq	%rax, _nze(%rip)
	movl	$d_one, %edx
	fld1
	fstpl	(%rdx)
	movl	$d_one+8, %eax
	fld1
	fchs
	fstpl	(%rax)
	fldz
	fchs
	fstpl	c_zero+8(%rip)
	fldz
	fstpl	c_zero(%rip)
	fldz
	fchs
	fstpl	c_zero+24(%rip)
	fldz
	fstpl	c_zero+16(%rip)
	movq	(%rdx), %rdx
	movq	%rdx, c_one(%rip)
	movq	(%rax), %rax
	movq	%rax, c_one+8(%rip)
	fldz
	fchs
	fstpl	c_one+24(%rip)
	fldz
	fstpl	c_one+16(%rip)
	flds	.LC8(%rip)
	fstpl	c_half+8(%rip)
	flds	.LC3(%rip)
	fstpl	c_half(%rip)
	fldz
	fchs
	fstpl	c_half+24(%rip)
	fldz
	fstpl	c_half+16(%rip)
	movsd	.LC9(%rip), %xmm1
	movsd	%xmm1, delta_int_double(%rip)
	movsd	.LC10(%rip), %xmm0
	movsd	%xmm0, delta_int_double+8(%rip)
	movsd	%xmm0, delta_int_double_neg(%rip)
	movsd	%xmm1, delta_int_double_neg+8(%rip)
	movl	$delta_blow, %edi
	call	_ZN10int_doubleC1Edd
	leaq	64(%rsp), %rdi
	movsd	.LC11(%rip), %xmm1
	movsd	.LC12(%rip), %xmm0
	call	_ZN10int_doubleC1Edd
	leaq	48(%rsp), %rdi
	movsd	.LC11(%rip), %xmm1
	movsd	.LC12(%rip), %xmm0
	call	_ZN10int_doubleC1Edd
	movq	48(%rsp), %rax
	movq	%rax, ln_gamma_err(%rip)
	movq	56(%rsp), %rax
	movq	%rax, ln_gamma_err+8(%rip)
	movl	$ln_gamma_err+16, %eax
	movq	64(%rsp), %rdx
	movq	%rdx, (%rax)
	movq	72(%rsp), %rdx
	movq	%rdx, 8(%rax)
	leaq	32(%rsp), %rdi
	movsd	.LC13(%rip), %xmm1
	movsd	.LC14(%rip), %xmm0
	call	_ZN10int_doubleC1Edd
	leaq	16(%rsp), %rdi
	movsd	.LC13(%rip), %xmm1
	movsd	.LC14(%rip), %xmm0
	call	_ZN10int_doubleC1Edd
	movq	16(%rsp), %rax
	movq	%rax, ln_gamma_err_2(%rip)
	movq	24(%rsp), %rax
	movq	%rax, ln_gamma_err_2+8(%rip)
	movl	$ln_gamma_err_2+16, %eax
	movq	32(%rsp), %rdx
	movq	%rdx, (%rax)
	movq	40(%rsp), %rdx
	movq	%rdx, 8(%rax)
	fldz
	fstpl	re_ln_gamma_err(%rip)
	fldz
	fchs
	fstpl	re_ln_gamma_err+8(%rip)
	addq	$88, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE1191:
	.size	_GLOBAL__I__Z16print_int_doubleRK10int_double, .-_GLOBAL__I__Z16print_int_doubleRK10int_double
	.section	.ctors,"aw",@progbits
	.align 8
	.quad	_GLOBAL__I__Z16print_int_doubleRK10int_double
	.section	.rodata.str1.8
	.align 8
.LC15:
	.string	"Error reading 64 bit part of zero from file. Exiting."
	.align 8
.LC16:
	.string	"Error reading 32 bit part of zero from file. Exiting."
	.align 8
.LC17:
	.string	"Error reading 8 bit part of zero from file. Exiting."
	.text
.globl _Z8get_zeroP8_IO_FILE
	.type	_Z8get_zeroP8_IO_FILE, @function
_Z8get_zeroP8_IO_FILE:
.LFB1186:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	pushq	%rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	subq	$24, %rsp
	.cfi_def_cfa_offset 48
	movq	%rdi, %rbx
	movq	%rsi, %rbp
	movq	%rsi, %rcx
	movl	$1, %edx
	movl	$8, %esi
	call	fread
	cmpq	$1, %rax
	je	.L14
	movl	$.LC15, %edi
	call	puts
	movl	$0, %edi
	call	exit
.L14:
	leaq	8(%rbx), %rdi
	movq	%rbp, %rcx
	movl	$1, %edx
	movl	$4, %esi
	call	fread
	cmpq	$1, %rax
	je	.L15
	movl	$.LC16, %edi
	call	puts
	movl	$0, %edi
	call	exit
.L15:
	leaq	15(%rsp), %rdi
	movq	%rbp, %rcx
	movl	$1, %edx
	movl	$1, %esi
	call	fread
	cmpq	$1, %rax
	je	.L16
	movl	$.LC17, %edi
	call	puts
	movl	$0, %edi
	call	exit
.L16:
	movzbl	15(%rsp), %eax
	movq	%rax, 16(%rbx)
	movq	%rbx, %rax
	addq	$24, %rsp
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE1186:
	.size	_Z8get_zeroP8_IO_FILE, .-_Z8get_zeroP8_IO_FILE
.globl _Z9_fpu_rnddv
	.type	_Z9_fpu_rnddv, @function
_Z9_fpu_rnddv:
.LFB1160:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	subq	$40, %rsp
	.cfi_def_cfa_offset 48
	call	crlibm_init
#APP
# 1497 "../includes/int_double12.0.h" 1
	stmxcsr old__SSE_cw(%rip)
# 0 "" 2
#NO_APP
	movl	old__SSE_cw(%rip), %eax
	andl	$8127, %eax
	orb	$32, %ah
	movl	%eax, new__SSE_cw(%rip)
#APP
# 1500 "../includes/int_double12.0.h" 1
	ldmxcsr new__SSE_cw(%rip)
# 0 "" 2
#NO_APP
	fldz
	fstl	d_zero_zero(%rip)
	fstl	d_zero_zero+8(%rip)
	fstl	d_zero(%rip)
	movq	_nze(%rip), %rax
	movq	%rax, d_zero+8(%rip)
	movq	%rax, d_neg_zero(%rip)
	movq	%rax, d_neg_zero+8(%rip)
	movq	%rax, d_neg_neg_zero(%rip)
	fstpl	d_neg_neg_zero+8(%rip)
	fldl	d_pi+8(%rip)
	fchs
	fstpl	24(%rsp)
	movsd	24(%rsp), %xmm0
	call	log_ru
	movsd	%xmm0, 24(%rsp)
	movsd	d_pi(%rip), %xmm0
	call	log_rd
	movsd	%xmm0, 8(%rsp)
	fldl	8(%rsp)
	fldl	24(%rsp)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	jbe	.L26
	movapd	%xmm0, %xmm1
	movsd	24(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L26:
	fldl	24(%rsp)
	fchs
	fstpl	d_ln_pi+8(%rip)
	movsd	%xmm0, d_ln_pi(%rip)
	fldl	d_two_pi+8(%rip)
	fchs
	fstpl	24(%rsp)
	movsd	24(%rsp), %xmm0
	call	log_ru
	movsd	%xmm0, 24(%rsp)
	movsd	d_two_pi(%rip), %xmm0
	call	log_rd
	movsd	%xmm0, 8(%rsp)
	fldl	8(%rsp)
	fldl	24(%rsp)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	jbe	.L27
	movapd	%xmm0, %xmm1
	movsd	24(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L27:
	fldl	24(%rsp)
	fchs
	fstpl	d_ln_two_pi+8(%rip)
	movsd	%xmm0, d_ln_two_pi(%rip)
	addq	$40, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE1160:
	.size	_Z9_fpu_rnddv, .-_Z9_fpu_rnddv
.globl _Z15set_h_bernoulliv
	.type	_Z15set_h_bernoulliv, @function
_Z15set_h_bernoulliv:
.LFB1161:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	subq	$64, %rsp
	.cfi_def_cfa_offset 72
	flds	.LC18(%rip)
	fstl	-96(%rsp)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup -96(%rsp),%xmm0
	movapd d_one(%rip),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,-120(%rsp)
# 0 "" 2
#NO_APP
	movq	-112(%rsp), %rax
	movq	%rax, bernoulli+8(%rip)
	movq	-120(%rsp), %rax
	movq	%rax, bernoulli(%rip)
	flds	.LC19(%rip)
	fstpl	-88(%rsp)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup -88(%rsp),%xmm0
	movapd bernoulli(%rip),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,-120(%rsp)
# 0 "" 2
#NO_APP
	movq	-112(%rsp), %rax
	movq	%rax, h_bernoulli+8(%rip)
	movq	-120(%rsp), %rax
	movq	%rax, h_bernoulli(%rip)
	flds	.LC20(%rip)
	fstpl	-104(%rsp)
#APP
# 653 "../includes/int_double12.0.h" 1
	movddup -104(%rsp),%xmm0
	movapd d_one(%rip),%xmm1
	shufpd $1,%xmm1,%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,-120(%rsp)
# 0 "" 2
#NO_APP
	movq	-112(%rsp), %rax
	movq	%rax, bernoulli+24(%rip)
	movl	$bernoulli+16, %eax
	movq	-120(%rsp), %rdx
	movq	%rdx, (%rax)
	flds	.LC21(%rip)
	fstpl	-80(%rsp)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup -80(%rsp),%xmm0
	movapd (%rax),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,-120(%rsp)
# 0 "" 2
#NO_APP
	movq	-112(%rsp), %rdx
	movq	%rdx, h_bernoulli+24(%rip)
	movq	-120(%rsp), %rdx
	movq	%rdx, h_bernoulli+16(%rip)
	flds	.LC22(%rip)
	fstpl	-72(%rsp)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup -72(%rsp),%xmm0
	movapd d_one(%rip),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,-120(%rsp)
# 0 "" 2
#NO_APP
	movq	-112(%rsp), %rdx
	movq	%rdx, bernoulli+40(%rip)
	movl	$bernoulli+32, %edx
	movq	-120(%rsp), %rcx
	movq	%rcx, (%rdx)
	flds	.LC23(%rip)
	fstpl	-64(%rsp)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup -64(%rsp),%xmm0
	movapd (%rdx),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,-120(%rsp)
# 0 "" 2
#NO_APP
	movq	-112(%rsp), %rdx
	movq	%rdx, h_bernoulli+40(%rip)
	movq	-120(%rsp), %rdx
	movq	%rdx, h_bernoulli+32(%rip)
	movl	$bernoulli+48, %edx
	movq	(%rax), %rcx
	movq	%rcx, (%rdx)
	movq	8(%rax), %rax
	movq	%rax, 8(%rdx)
	flds	.LC24(%rip)
	fstpl	-56(%rsp)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup -56(%rsp),%xmm0
	movapd (%rdx),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,-120(%rsp)
# 0 "" 2
#NO_APP
	movq	-112(%rsp), %rax
	movq	%rax, h_bernoulli+56(%rip)
	movq	-120(%rsp), %rax
	movq	%rax, h_bernoulli+48(%rip)
	flds	.LC25(%rip)
	fstpl	-48(%rsp)
	flds	.LC26(%rip)
	fstpl	-40(%rsp)
	flds	.LC27(%rip)
	fstpl	-32(%rsp)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup -48(%rsp),%xmm0
	movapd -40(%rsp),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,-120(%rsp)
# 0 "" 2
#NO_APP
	movq	-112(%rsp), %rax
	movq	%rax, bernoulli+72(%rip)
	movl	$bernoulli+64, %eax
	movq	-120(%rsp), %rdx
	movq	%rdx, (%rax)
	flds	.LC28(%rip)
	fstpl	-24(%rsp)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup -24(%rsp),%xmm0
	movapd (%rax),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,-120(%rsp)
# 0 "" 2
#NO_APP
	movq	-112(%rsp), %rax
	movq	%rax, h_bernoulli+72(%rip)
	movq	-120(%rsp), %rax
	movq	%rax, h_bernoulli+64(%rip)
	flds	.LC29(%rip)
	fstpl	-16(%rsp)
	flds	.LC30(%rip)
	fstpl	-8(%rsp)
	flds	.LC31(%rip)
	fstpl	(%rsp)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup -16(%rsp),%xmm0
	movapd -8(%rsp),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,-120(%rsp)
# 0 "" 2
#NO_APP
	movq	-112(%rsp), %rax
	movq	%rax, bernoulli+88(%rip)
	movl	$bernoulli+80, %eax
	movq	-120(%rsp), %rdx
	movq	%rdx, (%rax)
	flds	.LC32(%rip)
	fstl	8(%rsp)
	fxch	%st(1)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup 8(%rsp),%xmm0
	movapd (%rax),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,-120(%rsp)
# 0 "" 2
#NO_APP
	movq	-112(%rsp), %rax
	movq	%rax, h_bernoulli+88(%rip)
	movq	-120(%rsp), %rax
	movq	%rax, h_bernoulli+80(%rip)
	fstpl	16(%rsp)
	flds	.LC33(%rip)
	fstpl	24(%rsp)
	flds	.LC34(%rip)
	fstpl	32(%rsp)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup 16(%rsp),%xmm0
	movapd 24(%rsp),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,-120(%rsp)
# 0 "" 2
#NO_APP
	movq	-112(%rsp), %rax
	movq	%rax, bernoulli+104(%rip)
	movl	$bernoulli+96, %eax
	movq	-120(%rsp), %rdx
	movq	%rdx, (%rax)
	fstpl	40(%rsp)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup 40(%rsp),%xmm0
	movapd (%rax),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,-120(%rsp)
# 0 "" 2
#NO_APP
	movl	$h_bernoulli+104, %edx
	movq	-112(%rsp), %rax
	movq	%rax, (%rdx)
	movl	$h_bernoulli+96, %eax
	movq	-120(%rsp), %rcx
	movq	%rcx, (%rax)
	flds	.LC35(%rip)
	fstpl	48(%rsp)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup 48(%rsp),%xmm0
	movapd (%rax),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,-120(%rsp)
# 0 "" 2
#NO_APP
	movq	-112(%rsp), %rcx
	movq	%rcx, (%rdx)
	movq	-120(%rsp), %rdx
	movq	%rdx, (%rax)
	addq	$64, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE1161:
	.size	_Z15set_h_bernoulliv, .-_Z15set_h_bernoulliv
	.section	.rodata.str1.8
	.align 8
.LC37:
	.string	"Division by interval containing zero. Exiting"
	.text
.globl _Z8hurwitz0RK10int_double
	.type	_Z8hurwitz0RK10int_double, @function
_Z8hurwitz0RK10int_double:
.LFB1162:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	pushq	%rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	subq	$712, %rsp
	.cfi_def_cfa_offset 736
	movq	%rdi, %rbp
	cmpb	$0, h_bernoulli_initialised(%rip)
	jne	.L31
	call	_Z15set_h_bernoulliv
	movb	$1, h_bernoulli_initialised(%rip)
.L31:
	movq	$0, 144(%rsp)
	movq	$6, 152(%rsp)
	movq	$100, 160(%rsp)
	movq	$3528, 168(%rsp)
	movq	$219168, 176(%rsp)
	movq	$21257200, 184(%rsp)
	movl	$-1322081536, 192(%rsp)
	movl	$0, 196(%rsp)
	movq	$2, 80(%rsp)
	movq	$4, 88(%rsp)
	movq	$48, 96(%rsp)
	movq	$1440, 104(%rsp)
	movq	$80640, 112(%rsp)
	movq	$7257600, 120(%rsp)
	movq	$958003200, 128(%rsp)
	fldz
	fstpl	304(%rsp)
	fldz
	fchs
	fstpl	312(%rsp)
	movl	$0, %ebx
.L34:
	movl	%ebx, 44(%rsp)
	fildl	44(%rsp)
	fstpl	368(%rsp)
#APP
# 405 "../includes/int_double12.0.h" 1
	movddup 368(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd 0(%rbp),%XMM0
	movapd %XMM0,208(%rsp)
	
# 0 "" 2
#NO_APP
	fldl	208(%rsp)
	fstpl	32(%rsp)
	fldl	216(%rsp)
	fchs
	fstpl	24(%rsp)
	movsd	24(%rsp), %xmm0
	call	log_ru
	movsd	%xmm0, 24(%rsp)
	movsd	32(%rsp), %xmm0
	call	log_rd
	movsd	%xmm0, 8(%rsp)
	fldl	8(%rsp)
	fldl	24(%rsp)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	jbe	.L50
	movapd	%xmm0, %xmm1
	movsd	24(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L50:
	fldl	24(%rsp)
	fchs
	fstpl	360(%rsp)
	movsd	%xmm0, 352(%rsp)
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd 352(%rsp),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,208(%rsp)
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rax
	movq	%rax, 344(%rsp)
	movq	208(%rsp), %rax
	movq	%rax, 336(%rsp)
#APP
# 429 "../includes/int_double12.0.h" 1
	movapd 304(%rsp),%XMM0
	addpd 336(%rsp),%XMM0
	movapd %XMM0,304(%rsp)
	
# 0 "" 2
#NO_APP
	addl	$1, %ebx
	cmpl	$15, %ebx
	jne	.L34
	flds	.LC36(%rip)
	fstpl	376(%rsp)
#APP
# 405 "../includes/int_double12.0.h" 1
	movddup 376(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd 0(%rbp),%XMM0
	movapd %XMM0,208(%rsp)
	
# 0 "" 2
#NO_APP
	fldl	216(%rsp)
	fstl	296(%rsp)
	movq	208(%rsp), %rax
	movq	%rax, 288(%rsp)
	fchs
	fstpl	24(%rsp)
	movsd	24(%rsp), %xmm0
	call	log_ru
	movsd	%xmm0, 24(%rsp)
	movsd	288(%rsp), %xmm0
	call	log_rd
	movsd	%xmm0, 8(%rsp)
	fldl	8(%rsp)
	fldl	24(%rsp)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	jbe	.L51
	movapd	%xmm0, %xmm1
	movsd	24(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L51:
	fldl	24(%rsp)
	fchs
	fstpl	280(%rsp)
	movsd	%xmm0, 272(%rsp)
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd 272(%rsp),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,208(%rsp)
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rax
	movq	%rax, 264(%rsp)
	movq	208(%rsp), %rax
	movq	%rax, 256(%rsp)
	flds	.LC3(%rip)
	fstpl	424(%rsp)
#APP
# 579 "../includes/int_double12.0.h" 1
	movddup 424(%rsp),%xmm0
	movapd %xmm0,%xmm1
	xorpd d_neg_zero(%rip),%xmm1
	movapd 256(%rsp),%xmm2
	mulpd %xmm2,%xmm0
	mulpd %xmm2,%xmm1
	shufpd $1,%xmm1,%xmm1
	minpd %xmm1,%xmm0
	movapd %xmm0,208(%rsp)
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rax
	movq	%rax, 408(%rsp)
	movq	208(%rsp), %rax
	movq	%rax, 400(%rsp)
	flds	.LC19(%rip)
	fstl	472(%rsp)
#APP
# 579 "../includes/int_double12.0.h" 1
	movddup 472(%rsp),%xmm0
	movapd %xmm0,%xmm1
	xorpd d_neg_zero(%rip),%xmm1
	movapd 0(%rbp),%xmm2
	mulpd %xmm2,%xmm0
	mulpd %xmm2,%xmm1
	shufpd $1,%xmm1,%xmm1
	minpd %xmm1,%xmm0
	movapd %xmm0,208(%rsp)
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rax
	movq	%rax, 456(%rsp)
	movq	208(%rsp), %rax
	movq	%rax, 448(%rsp)
	flds	.LC20(%rip)
	fstpl	504(%rsp)
	fstpl	568(%rsp)
#APP
# 579 "../includes/int_double12.0.h" 1
	movddup 568(%rsp),%xmm0
	movapd %xmm0,%xmm1
	xorpd d_neg_zero(%rip),%xmm1
	movapd 288(%rsp),%xmm2
	mulpd %xmm2,%xmm0
	mulpd %xmm2,%xmm1
	shufpd $1,%xmm1,%xmm1
	minpd %xmm1,%xmm0
	movapd %xmm0,208(%rsp)
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rax
	movq	%rax, 552(%rsp)
	movq	208(%rsp), %rax
	movq	%rax, 544(%rsp)
#APP
# 524 "../includes/int_double12.0.h" 1
	movapd 272(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 544(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,208(%rsp)
	
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rax
	movq	%rax, 536(%rsp)
	movq	208(%rsp), %rax
	movq	%rax, 528(%rsp)
#APP
# 524 "../includes/int_double12.0.h" 1
	movapd 256(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 288(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,208(%rsp)
	
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rax
	movq	%rax, 600(%rsp)
	movq	208(%rsp), %rax
	movq	%rax, 592(%rsp)
#APP
# 454 "../includes/int_double12.0.h" 1
	movapd 592(%rsp),%xmm0
	shufpd $1,%xmm0,%xmm0
	addpd 256(%rsp),%xmm0
	movapd %xmm0,208(%rsp)
	
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rax
	movq	%rax, 584(%rsp)
	movq	208(%rsp), %rax
	movq	%rax, 576(%rsp)
#APP
# 392 "../includes/int_double12.0.h" 1
	movapd 576(%rsp),%XMM0
	addpd 528(%rsp),%XMM0
	movapd %XMM0,208(%rsp)
	
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rax
	movq	%rax, 520(%rsp)
	movq	208(%rsp), %rax
	movq	%rax, 512(%rsp)
#APP
# 467 "../includes/int_double12.0.h" 1
	movddup 504(%rsp),%XMM0
	movapd 512(%rsp),%XMM1
	addsubpd %XMM0,%XMM1
	movapd %XMM1,208(%rsp)
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rax
	movq	%rax, 488(%rsp)
	movq	208(%rsp), %rax
	movq	%rax, 480(%rsp)
#APP
# 454 "../includes/int_double12.0.h" 1
	movapd 448(%rsp),%xmm0
	shufpd $1,%xmm0,%xmm0
	addpd 480(%rsp),%xmm0
	movapd %xmm0,208(%rsp)
	
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rax
	movq	%rax, 440(%rsp)
	movq	208(%rsp), %rax
	movq	%rax, 432(%rsp)
#APP
# 454 "../includes/int_double12.0.h" 1
	movapd 400(%rsp),%xmm0
	shufpd $1,%xmm0,%xmm0
	addpd 432(%rsp),%xmm0
	movapd %xmm0,208(%rsp)
	
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rax
	movq	%rax, 392(%rsp)
	movq	208(%rsp), %rax
	movq	%rax, 384(%rsp)
#APP
# 429 "../includes/int_double12.0.h" 1
	movapd 304(%rsp),%XMM0
	addpd 384(%rsp),%XMM0
	movapd %XMM0,304(%rsp)
	
# 0 "" 2
# 1003 "../includes/int_double12.0.h" 1
	movapd 288(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,332(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 332(%rsp)
	je	.L37
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	288(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L37:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 288(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd d_one(%rip),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,208(%rsp)
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rax
	movq	%rax, 296(%rsp)
	movq	208(%rsp), %rax
	movq	%rax, 288(%rsp)
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd 288(%rsp),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,208(%rsp)
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rax
	movq	%rax, 248(%rsp)
	movq	208(%rsp), %rax
	movq	%rax, 240(%rsp)
	movl	$0, %eax
	movl	$0, %edx
	leaq	80(%rsp), %rbx
	flds	.LC38(%rip)
	leaq	144(%rsp), %rcx
	fld	%st(0)
	movl	$h_bernoulli, %edi
.L40:
	fildll	(%rbx,%rax)
	cmpq	$0, (%rbx,%rax)
	jns	.L38
	fadd	%st(2), %st
.L38:
	fstpl	664(%rsp)
#APP
# 579 "../includes/int_double12.0.h" 1
	movddup 664(%rsp),%xmm0
	movapd %xmm0,%xmm1
	xorpd d_neg_zero(%rip),%xmm1
	movapd 272(%rsp),%xmm2
	mulpd %xmm2,%xmm0
	mulpd %xmm2,%xmm1
	shufpd $1,%xmm1,%xmm1
	minpd %xmm1,%xmm0
	movapd %xmm0,208(%rsp)
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rsi
	movq	%rsi, 648(%rsp)
	movq	208(%rsp), %rsi
	movq	%rsi, 640(%rsp)
	fildll	(%rcx,%rax)
	cmpq	$0, (%rcx,%rax)
	jns	.L39
	fadd	%st(1), %st
.L39:
	fstpl	72(%rsp)
	fldl	72(%rsp)
	fstl	672(%rsp)
	fchs
	fstpl	680(%rsp)
#APP
# 454 "../includes/int_double12.0.h" 1
	movapd 640(%rsp),%xmm0
	shufpd $1,%xmm0,%xmm0
	addpd 672(%rsp),%xmm0
	movapd %xmm0,208(%rsp)
	
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rsi
	movq	%rsi, 632(%rsp)
	movq	208(%rsp), %rsi
	movq	%rsi, 624(%rsp)
	movslq	%edx, %rsi
	salq	$4, %rsi
#APP
# 524 "../includes/int_double12.0.h" 1
	movapd 624(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd (%rdi,%rsi),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,208(%rsp)
	
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rsi
	movq	%rsi, 616(%rsp)
	movq	208(%rsp), %rsi
	movq	%rsi, 608(%rsp)
#APP
# 524 "../includes/int_double12.0.h" 1
	movapd 288(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 608(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,208(%rsp)
	
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rsi
	movq	%rsi, 232(%rsp)
	movq	208(%rsp), %rsi
	movq	%rsi, 224(%rsp)
#APP
# 429 "../includes/int_double12.0.h" 1
	movapd 304(%rsp),%XMM0
	addpd 224(%rsp),%XMM0
	movapd %XMM0,304(%rsp)
	
# 0 "" 2
# 524 "../includes/int_double12.0.h" 1
	movapd 240(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 288(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,208(%rsp)
	
# 0 "" 2
#NO_APP
	movq	216(%rsp), %rsi
	movq	%rsi, 296(%rsp)
	movq	208(%rsp), %rsi
	movq	%rsi, 288(%rsp)
	addl	$1, %edx
	addq	$8, %rax
	cmpl	$5, %edx
	jne	.L40
	fstp	%st(0)
	fstp	%st(0)
	fldl	224(%rsp)
	fldl	232(%rsp)
	fucomi	%st(1), %st
	jb	.L52
	fstp	%st(0)
	fstpl	232(%rsp)
	jmp	.L43
.L52:
	fstp	%st(1)
	fstpl	224(%rsp)
.L43:
#APP
# 429 "../includes/int_double12.0.h" 1
	movapd 304(%rsp),%XMM0
	addpd 224(%rsp),%XMM0
	movapd %XMM0,304(%rsp)
	
# 0 "" 2
#NO_APP
	movq	304(%rsp), %rax
	movq	%rax, 688(%rsp)
	movq	312(%rsp), %rax
	movq	%rax, 696(%rsp)
	movsd	688(%rsp), %xmm0
	movsd	696(%rsp), %xmm1
	addq	$712, %rsp
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE1162:
	.size	_Z8hurwitz0RK10int_double, .-_Z8hurwitz0RK10int_double
	.section	.rodata.str1.8
	.align 8
.LC40:
	.string	"Endpoints more than Pi/2 apart in sin_cos. Exiting."
	.align 8
.LC41:
	.string	"Weird error in sin_cos, mask = %d Exiting.\n"
	.align 8
.LC42:
	.string	"\ncos_left: %20.18e\ncos_right: %20.18e\n"
	.align 8
.LC43:
	.string	"sin_left: %20.18e\nsin_right: %20.18e\n"
	.section	.text._Z3powRK10int_doubleRK11int_complex,"axG",@progbits,_Z3powRK10int_doubleRK11int_complex,comdat
	.weak	_Z3powRK10int_doubleRK11int_complex
	.type	_Z3powRK10int_doubleRK11int_complex, @function
_Z3powRK10int_doubleRK11int_complex:
.LFB1151:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	pushq	%r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	pushq	%rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	subq	$160, %rsp
	.cfi_def_cfa_offset 192
	movq	%rdi, %rbx
	movq	%rsi, %rbp
	movq	%rdx, %r12
	fldl	8(%rsi)
	fchs
	fstpl	16(%rsp)
	movsd	16(%rsp), %xmm0
	call	log_ru
	movsd	%xmm0, 16(%rsp)
	movsd	0(%rbp), %xmm0
	call	log_rd
	movsd	%xmm0, 8(%rsp)
	fldl	8(%rsp)
	fldl	16(%rsp)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	jbe	.L104
	movapd	%xmm0, %xmm1
	movsd	16(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L104:
	fldl	16(%rsp)
	fchs
	fstpl	104(%rsp)
	movsd	%xmm0, 96(%rsp)
#APP
# 524 "../includes/int_double12.0.h" 1
	movapd 16(%r12),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 96(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,80(%rsp)
	
# 0 "" 2
#NO_APP
	fldl	88(%rsp)
	fldl	80(%rsp)
	fxch	%st(1)
	fstl	136(%rsp)
	fxch	%st(1)
	fstl	128(%rsp)
	faddp	%st, %st(1)
	fldl	.LC39(%rip)
	fucomip	%st(1), %st
	fstp	%st(0)
	jb	.L105
	movl	$.LC40, %edi
	call	puts
	movl	$1, %edi
	call	exit
.L105:
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd d_pi(%rip),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,120(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 120(%rsp)
	je	.L58
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	movl	$d_pi, %edi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L58:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd d_pi(%rip),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 128(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,80(%rsp)
# 0 "" 2
#NO_APP
	fldl	88(%rsp)
	fstpl	16(%rsp)
	fldl	80(%rsp)
	fstpl	24(%rsp)
	movsd	24(%rsp), %xmm0
	call	cospi_rd
	movsd	%xmm0, 32(%rsp)
	fldl	16(%rsp)
	fchs
	fstpl	16(%rsp)
	movsd	16(%rsp), %xmm0
	call	cospi_rd
	movsd	%xmm0, 40(%rsp)
	movsd	24(%rsp), %xmm0
	call	sinpi_rd
	movsd	%xmm0, 24(%rsp)
	movsd	16(%rsp), %xmm0
	call	sinpi_rd
	movsd	%xmm0, 16(%rsp)
	fldz
	fldl	32(%rsp)
	fucomip	%st(1), %st
	sbbl	%esi, %esi
	notl	%esi
	andl	$8, %esi
	movl	%esi, %eax
	orl	$4, %eax
	fldl	40(%rsp)
	fucomip	%st(1), %st
	cmovae	%eax, %esi
	movl	%esi, %eax
	orl	$2, %eax
	fldl	24(%rsp)
	fucomip	%st(1), %st
	cmovae	%eax, %esi
	leal	1(%rsi), %eax
	fldl	16(%rsp)
	fucomip	%st(1), %st
	fstp	%st(0)
	cmovae	%eax, %esi
	cmpb	$15, %sil
	ja	.L68
	movzbl	%sil, %eax
	jmp	*.L77(,%rax,8)
	.section	.rodata._Z3powRK10int_doubleRK11int_complex,"aG",@progbits,_Z3powRK10int_doubleRK11int_complex,comdat
	.align 8
	.align 4
.L77:
	.quad	.L69
	.quad	.L68
	.quad	.L70
	.quad	.L71
	.quad	.L72
	.quad	.L68
	.quad	.L68
	.quad	.L68
	.quad	.L68
	.quad	.L68
	.quad	.L68
	.quad	.L73
	.quad	.L74
	.quad	.L75
	.quad	.L68
	.quad	.L76
	.section	.text._Z3powRK10int_doubleRK11int_complex,"axG",@progbits,_Z3powRK10int_doubleRK11int_complex,comdat
.L69:
	fldl	16(%rsp)
	fstpl	64(%rsp)
	fldl	24(%rsp)
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rax
	movq	%rax, 72(%rsp)
	fldl	32(%rsp)
	fstpl	48(%rsp)
	fldl	40(%rsp)
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rax
	movq	%rax, 56(%rsp)
	jmp	.L78
.L70:
	fldl	16(%rsp)
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double_neg(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	80(%rsp), %rax
	movq	%rax, 64(%rsp)
	fldl	24(%rsp)
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rax
	movq	%rax, 72(%rsp)
	fld1
	fchs
	fstpl	48(%rsp)
	fldl	32(%rsp)
	fldl	40(%rsp)
	fxch	%st(1)
	fucomi	%st(1), %st
	fstp	%st(1)
	jb	.L106
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rax
	movq	%rax, 56(%rsp)
	jmp	.L78
.L106:
	fstp	%st(0)
	fldl	40(%rsp)
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rax
	movq	%rax, 56(%rsp)
	jmp	.L78
.L71:
	fldl	16(%rsp)
	fstpl	64(%rsp)
	fldl	24(%rsp)
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rax
	movq	%rax, 72(%rsp)
	fldl	40(%rsp)
	fstpl	48(%rsp)
	fldl	32(%rsp)
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rax
	movq	%rax, 56(%rsp)
	jmp	.L78
.L72:
	fld1
	fchs
	fstpl	64(%rsp)
	fldl	16(%rsp)
	fldl	24(%rsp)
	fxch	%st(1)
	fucomi	%st(1), %st
	fstp	%st(1)
	jb	.L107
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rax
	movq	%rax, 72(%rsp)
	jmp	.L83
.L107:
	fstp	%st(0)
	fldl	24(%rsp)
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rax
	movq	%rax, 72(%rsp)
.L83:
	fldl	32(%rsp)
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double_neg(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	80(%rsp), %rax
	movq	%rax, 48(%rsp)
	fldl	40(%rsp)
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rax
	movq	%rax, 56(%rsp)
	jmp	.L78
.L73:
	fldl	16(%rsp)
	fldl	24(%rsp)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	movq	16(%rsp), %rax
	cmovae	24(%rsp), %rax
	movq	%rax, 64(%rsp)
	fld1
	fchs
	fstpl	72(%rsp)
	fldl	40(%rsp)
	fstpl	48(%rsp)
	fldl	32(%rsp)
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rax
	movq	%rax, 56(%rsp)
	jmp	.L78
.L74:
	fldl	24(%rsp)
	fstpl	64(%rsp)
	fldl	16(%rsp)
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rax
	movq	%rax, 72(%rsp)
	fldl	32(%rsp)
	fstpl	48(%rsp)
	fldl	40(%rsp)
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rax
	movq	%rax, 56(%rsp)
	jmp	.L78
.L75:
	fld1
	fchs
	fstpl	56(%rsp)
	fldl	40(%rsp)
	fldl	32(%rsp)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	movq	40(%rsp), %rax
	cmovae	32(%rsp), %rax
	movq	%rax, 48(%rsp)
	fldl	24(%rsp)
	fstpl	64(%rsp)
	fldl	16(%rsp)
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rax
	movq	%rax, 72(%rsp)
	jmp	.L78
.L76:
	fldl	24(%rsp)
	fstpl	64(%rsp)
	fldl	16(%rsp)
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rax
	movq	%rax, 72(%rsp)
	fldl	40(%rsp)
	fstpl	48(%rsp)
	fldl	32(%rsp)
	fstpl	120(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 120(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rax
	movq	%rax, 56(%rsp)
	jmp	.L78
.L68:
	movzbl	%sil, %esi
	movl	$.LC41, %edi
	movl	$0, %eax
	call	printf
	movsd	40(%rsp), %xmm1
	movsd	32(%rsp), %xmm0
	movl	$.LC42, %edi
	movl	$2, %eax
	call	printf
	movsd	16(%rsp), %xmm1
	movsd	24(%rsp), %xmm0
	movl	$.LC43, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L78:
#APP
# 524 "../includes/int_double12.0.h" 1
	movapd (%r12),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 96(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,80(%rsp)
	
# 0 "" 2
#NO_APP
	fldl	80(%rsp)
	fstpl	32(%rsp)
	fldl	88(%rsp)
	fchs
	fstpl	16(%rsp)
	movsd	16(%rsp), %xmm0
	call	exp_ru
	movsd	%xmm0, 16(%rsp)
	movsd	32(%rsp), %xmm0
	call	exp_rd
	movsd	%xmm0, 8(%rsp)
	fldl	8(%rsp)
	fldl	16(%rsp)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	jbe	.L108
	movapd	%xmm0, %xmm1
	movsd	16(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L108:
	fldl	16(%rsp)
	fchs
	fstpl	152(%rsp)
	movsd	%xmm0, 144(%rsp)
#APP
# 524 "../includes/int_double12.0.h" 1
	movapd 144(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 64(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rdx
	movq	80(%rsp), %rax
#APP
# 524 "../includes/int_double12.0.h" 1
	movapd 144(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 48(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,80(%rsp)
	
# 0 "" 2
#NO_APP
	movq	88(%rsp), %rcx
	movq	%rcx, 8(%rbx)
	movq	80(%rsp), %rcx
	movq	%rcx, (%rbx)
	movq	%rdx, 24(%rbx)
	movq	%rax, 16(%rbx)
	movq	%rbx, %rax
	addq	$160, %rsp
	.cfi_def_cfa_offset 32
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE1151:
	.size	_Z3powRK10int_doubleRK11int_complex, .-_Z3powRK10int_doubleRK11int_complex
	.section	.rodata.str1.8
	.align 8
.LC44:
	.string	"Error in atan2. Can't handle arguments on negative real axis. Exiting."
	.align 8
.LC45:
	.string	"Error in atan2. Both x and y contain zero. Exiting."
	.align 8
.LC46:
	.string	"Bad mask in atan2. Mask=%d. Exiting.\n"
	.section	.rodata.str1.1
.LC47:
	.string	"x="
.LC48:
	.string	"y="
	.section	.rodata.str1.8
	.align 8
.LC49:
	.string	"Division by zero in int_double / double. Exiting."
	.section	.text._Z3logRK11int_complex,"axG",@progbits,_Z3logRK11int_complex,comdat
	.weak	_Z3logRK11int_complex
	.type	_Z3logRK11int_complex, @function
_Z3logRK11int_complex:
.LFB1157:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	pushq	%r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	pushq	%rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	subq	$192, %rsp
	.cfi_def_cfa_offset 224
	movq	%rdi, %rbp
	movq	%rsi, %rbx
	leaq	16(%rsi), %r12
	fldl	24(%rsi)
	fldl	16(%rsi)
	fldl	8(%rsi)
	fldl	(%rsi)
	fldz
	fucomi	%st(1), %st
	fstp	%st(1)
	fxch	%st(1)
	seta	%sil
	movzbl	%sil, %esi
	sall	$3, %esi
	fchs
	fxch	%st(1)
	movl	%esi, %eax
	orl	$4, %eax
	fucomi	%st(1), %st
	fstp	%st(1)
	cmova	%eax, %esi
	movl	%esi, %eax
	orl	$2, %eax
	fucomi	%st(1), %st
	fstp	%st(1)
	fxch	%st(1)
	cmova	%eax, %esi
	fchs
	fxch	%st(1)
	movl	%esi, %eax
	orl	$1, %eax
	fucomip	%st(1), %st
	fstp	%st(0)
	cmova	%eax, %esi
	cmpl	$15, %esi
	ja	.L119
	mov	%esi, %eax
	jmp	*.L125(,%rax,8)
	.section	.rodata._Z3logRK11int_complex,"aG",@progbits,_Z3logRK11int_complex,comdat
	.align 8
	.align 4
.L125:
	.quad	.L120
	.quad	.L119
	.quad	.L120
	.quad	.L120
	.quad	.L119
	.quad	.L119
	.quad	.L119
	.quad	.L119
	.quad	.L121
	.quad	.L119
	.quad	.L122
	.quad	.L123
	.quad	.L121
	.quad	.L119
	.quad	.L124
	.quad	.L123
	.section	.text._Z3logRK11int_complex,"axG",@progbits,_Z3logRK11int_complex,comdat
.L120:
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd (%rbx),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,64(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 64(%rsp)
	je	.L126
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	movq	%rbx, %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L126:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd (%rbx),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd (%r12),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,48(%rsp)
# 0 "" 2
#NO_APP
	fldl	48(%rsp)
	fstpl	16(%rsp)
	fldl	56(%rsp)
	fchs
	fstpl	24(%rsp)
	movsd	24(%rsp), %xmm0
	call	atan_ru
	movsd	%xmm0, 24(%rsp)
	movsd	16(%rsp), %xmm0
	call	atan_rd
	movsd	%xmm0, 32(%rsp)
	fldl	32(%rsp)
	fldl	24(%rsp)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	jbe	.L154
	movapd	%xmm0, %xmm1
	movsd	24(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L154:
	fldl	24(%rsp)
	fchs
	fstpl	24(%rsp)
	jmp	.L129
.L121:
#APP
# 479 "../includes/int_double12.0.h" 1
	movapd (%rbx),%XMM0
	shufpd $1,%XMM0,%XMM0
	movupd %XMM0,48(%rsp)
	
# 0 "" 2
#NO_APP
	movq	56(%rsp), %rax
	movq	%rax, 120(%rsp)
	movq	48(%rsp), %rax
	movq	%rax, 112(%rsp)
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd (%r12),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,64(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 64(%rsp)
	je	.L130
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	movq	%r12, %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L130:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd (%r12),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 112(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,48(%rsp)
# 0 "" 2
#NO_APP
	fldl	48(%rsp)
	fstpl	24(%rsp)
	fldl	56(%rsp)
	fchs
	fstpl	16(%rsp)
	movsd	16(%rsp), %xmm0
	call	atan_ru
	movsd	%xmm0, 16(%rsp)
	movsd	24(%rsp), %xmm0
	call	atan_rd
	movsd	%xmm0, 8(%rsp)
	fldl	8(%rsp)
	fldl	16(%rsp)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	jbe	.L155
	movapd	%xmm0, %xmm1
	movsd	16(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L155:
	fldl	16(%rsp)
	fchs
	fstpl	104(%rsp)
	movsd	%xmm0, 96(%rsp)
#APP
# 392 "../includes/int_double12.0.h" 1
	movapd 96(%rsp),%XMM0
	addpd d_pi_2(%rip),%XMM0
	movapd %XMM0,48(%rsp)
	
# 0 "" 2
#NO_APP
	fldl	56(%rsp)
	fstpl	24(%rsp)
	fldl	48(%rsp)
	fstpl	32(%rsp)
	jmp	.L129
.L123:
#APP
# 479 "../includes/int_double12.0.h" 1
	movapd (%r12),%XMM0
	shufpd $1,%XMM0,%XMM0
	movupd %XMM0,48(%rsp)
	
# 0 "" 2
#NO_APP
	movq	56(%rsp), %rax
	movq	%rax, 152(%rsp)
	movq	48(%rsp), %rax
	movq	%rax, 144(%rsp)
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd 144(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,64(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 64(%rsp)
	je	.L133
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	144(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L133:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 144(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd (%rbx),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,48(%rsp)
# 0 "" 2
#NO_APP
	fldl	48(%rsp)
	fstpl	24(%rsp)
	fldl	56(%rsp)
	fchs
	fstpl	16(%rsp)
	movsd	16(%rsp), %xmm0
	call	atan_ru
	movsd	%xmm0, 16(%rsp)
	movsd	24(%rsp), %xmm0
	call	atan_rd
	movsd	%xmm0, 8(%rsp)
	fldl	8(%rsp)
	fldl	16(%rsp)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	jbe	.L156
	movapd	%xmm0, %xmm1
	movsd	16(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L156:
	fldl	16(%rsp)
	fchs
	fstpl	136(%rsp)
	movsd	%xmm0, 128(%rsp)
#APP
# 454 "../includes/int_double12.0.h" 1
	movapd d_pi_2(%rip),%xmm0
	shufpd $1,%xmm0,%xmm0
	addpd 128(%rsp),%xmm0
	movapd %xmm0,48(%rsp)
	
# 0 "" 2
#NO_APP
	fldl	56(%rsp)
	fstpl	24(%rsp)
	fldl	48(%rsp)
	fstpl	32(%rsp)
	jmp	.L129
.L124:
	movl	$.LC44, %edi
	call	puts
	movl	$1, %edi
	call	exit
.L122:
	movl	$.LC45, %edi
	call	puts
	movl	$1, %edi
	call	exit
.L119:
	movl	$.LC46, %edi
	movl	$0, %eax
	call	printf
	movl	$.LC47, %edi
	movl	$0, %eax
	call	printf
	movl	$32, %edi
	call	putchar
	movq	%rbx, %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$.LC48, %edi
	movl	$0, %eax
	call	printf
	movl	$32, %edi
	call	putchar
	movq	%r12, %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L129:
	flds	.LC19(%rip)
	fstpl	72(%rsp)
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd (%r12),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,48(%rsp)
# 0 "" 2
#NO_APP
	movq	56(%rsp), %rax
	movq	%rax, 184(%rsp)
	movq	48(%rsp), %rax
	movq	%rax, 176(%rsp)
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd (%rbx),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,48(%rsp)
# 0 "" 2
#NO_APP
	movq	56(%rsp), %rax
	movq	%rax, 168(%rsp)
	movq	48(%rsp), %rax
	movq	%rax, 160(%rsp)
#APP
# 392 "../includes/int_double12.0.h" 1
	movapd 160(%rsp),%XMM0
	addpd 176(%rsp),%XMM0
	movapd %XMM0,48(%rsp)
	
# 0 "" 2
#NO_APP
	fldl	48(%rsp)
	fstpl	40(%rsp)
	fldl	56(%rsp)
	fchs
	fstpl	16(%rsp)
	movsd	16(%rsp), %xmm0
	call	log_ru
	movsd	%xmm0, 16(%rsp)
	movsd	40(%rsp), %xmm0
	call	log_rd
	movsd	%xmm0, 8(%rsp)
	fldl	8(%rsp)
	fldl	16(%rsp)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	jbe	.L157
	movapd	%xmm0, %xmm1
	movsd	16(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L157:
	fldl	16(%rsp)
	fchs
	fstpl	88(%rsp)
	movsd	%xmm0, 80(%rsp)
	fldl	72(%rsp)
	fldz
	fxch	%st(1)
	fucomi	%st(1), %st
	fstp	%st(1)
	jbe	.L158
	fstp	%st(0)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup 72(%rsp),%xmm0
	movapd 80(%rsp),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,48(%rsp)
# 0 "" 2
#NO_APP
	movq	56(%rsp), %rdx
	movq	48(%rsp), %rax
	jmp	.L140
.L158:
	fldz
	fucomip	%st(1), %st
	jbe	.L159
	fchs
	fstpl	64(%rsp)
#APP
# 653 "../includes/int_double12.0.h" 1
	movddup 64(%rsp),%xmm0
	movapd 80(%rsp),%xmm1
	shufpd $1,%xmm1,%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,48(%rsp)
# 0 "" 2
#NO_APP
	movq	56(%rsp), %rdx
	movq	48(%rsp), %rax
	jmp	.L140
.L159:
	fstp	%st(0)
	movl	$.LC49, %edi
	call	puts
	movl	$1, %edi
	call	exit
.L140:
	movq	%rdx, 8(%rbp)
	movq	%rax, 0(%rbp)
	fldl	24(%rsp)
	fstpl	24(%rbp)
	fldl	32(%rsp)
	fstpl	16(%rbp)
	movq	%rbp, %rax
	addq	$192, %rsp
	.cfi_def_cfa_offset 32
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE1157:
	.size	_Z3logRK11int_complex, .-_Z3logRK11int_complex
	.text
.globl _Z1S10int_doublem
	.type	_Z1S10int_doublem, @function
_Z1S10int_doublem:
.LFB1187:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	pushq	%r13
	.cfi_def_cfa_offset 16
	.cfi_offset 13, -16
	pushq	%r12
	.cfi_def_cfa_offset 24
	.cfi_offset 12, -24
	pushq	%rbp
	.cfi_def_cfa_offset 32
	.cfi_offset 6, -32
	pushq	%rbx
	.cfi_def_cfa_offset 40
	.cfi_offset 3, -40
	subq	$2120, %rsp
	.cfi_def_cfa_offset 2160
	movq	%rdi, %r13
	movsd	%xmm0, 48(%rsp)
	movsd	%xmm1, 56(%rsp)
	flds	.LC3(%rip)
	fstpl	296(%rsp)
#APP
# 561 "../includes/int_double12.0.h" 1
	movddup 296(%rsp),%xmm0
	movapd %xmm0,%xmm1
	xorpd d_neg_zero(%rip),%xmm1
	movapd 48(%rsp),%xmm2
	mulpd %xmm2,%xmm0
	mulpd %xmm2,%xmm1
	shufpd $1,%xmm1,%xmm1
	minpd %xmm1,%xmm0
	movapd %xmm0,256(%rsp)
# 0 "" 2
#NO_APP
	flds	.LC50(%rip)
	fstpl	232(%rsp)
	flds	.LC51(%rip)
	fstpl	224(%rsp)
	movq	264(%rsp), %rax
	movq	%rax, 248(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 240(%rsp)
	movq	224(%rsp), %rax
	movq	%rax, 64(%rsp)
	movq	232(%rsp), %rax
	movq	%rax, 72(%rsp)
	movq	240(%rsp), %rax
	movq	%rax, 80(%rsp)
	movq	248(%rsp), %rax
	movq	%rax, 88(%rsp)
	cmpb	$0, h_bernoulli_initialised(%rip)
	jne	.L161
	call	_Z15set_h_bernoulliv
	movb	$1, h_bernoulli_initialised(%rip)
.L161:
	fldl	80(%rsp)
	flds	.LC52(%rip)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	jb	.L191
	flds	.LC3(%rip)
	fstpl	1640(%rsp)
#APP
# 561 "../includes/int_double12.0.h" 1
	movddup 1640(%rsp),%xmm0
	movapd %xmm0,%xmm1
	xorpd d_neg_zero(%rip),%xmm1
	movapd d_ln_two_pi(%rip),%xmm2
	mulpd %xmm2,%xmm0
	mulpd %xmm2,%xmm1
	shufpd $1,%xmm1,%xmm1
	minpd %xmm1,%xmm0
	movapd %xmm0,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 1656(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 1648(%rsp)
	leaq	1536(%rsp), %rdi
	leaq	64(%rsp), %rsi
	call	_Z3logRK11int_complex
	flds	.LC3(%rip)
	fstpl	1496(%rsp)
#APP
# 1187 "../includes/int_double12.0.h" 1
	movddup 1496(%rsp),%xmm0
	movapd 64(%rsp),%xmm1
	addsubpd %xmm0,%xmm1
	movapd %xmm1,1504(%rsp)
	movapd 80(%rsp),%xmm1
	movapd %xmm1,1520(%rsp)
	
# 0 "" 2
# 1318 "../includes/int_double12.0.h" 1
	movapd 1520(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1552(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	shufpd $1,%xmm7,%xmm7
	movapd 1504(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1536(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm7,%xmm0
	movapd %xmm0,1568(%rsp)
	movapd 1520(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1536(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	movapd 1504(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1552(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm0,%xmm7
	movapd %xmm7,1584(%rsp)
# 0 "" 2
# 1157 "../includes/int_double12.0.h" 1
	movapd 64(%rsp),%xmm0
	shufpd $1,%xmm0,%xmm0
	addpd 1568(%rsp),%xmm0
	movapd %xmm0,1600(%rsp)
	movapd 80(%rsp),%xmm1
	shufpd $1,%xmm1,%xmm1
	addpd 1584(%rsp),%xmm1
	movapd %xmm1,1616(%rsp)
	
# 0 "" 2
# 1062 "../includes/int_double12.0.h" 1
	movapd 1600(%rsp),%xmm0
	addpd 1648(%rsp),%xmm0
	movapd %xmm0,656(%rsp)
	movapd 1616(%rsp),%xmm1
	movapd %xmm1,672(%rsp)
	
# 0 "" 2
#NO_APP
	movq	656(%rsp), %rax
	movq	%rax, 160(%rsp)
	movq	664(%rsp), %rax
	movq	%rax, 168(%rsp)
	movq	672(%rsp), %rax
	movq	%rax, 176(%rsp)
	movq	680(%rsp), %rax
	movq	%rax, 184(%rsp)
	movq	bernoulli(%rip), %rax
	movq	%rax, 1376(%rsp)
	movq	bernoulli+8(%rip), %rax
	movq	%rax, 1384(%rsp)
	movq	d_zero(%rip), %rax
	movq	%rax, 1392(%rsp)
	movq	d_zero+8(%rip), %rax
	movq	%rax, 1400(%rsp)
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd 80(%rsp),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 1672(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 1664(%rsp)
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd 64(%rsp),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 1688(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 1680(%rsp)
#APP
# 392 "../includes/int_double12.0.h" 1
	movapd 1680(%rsp),%XMM0
	addpd 1664(%rsp),%XMM0
	movapd %XMM0,256(%rsp)
	
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 1768(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 1760(%rsp)
#APP
# 1033 "../includes/int_double12.0.h" 1
	movapd 64(%rsp),%xmm0
	movapd %xmm0,1696(%rsp)
	movapd 80(%rsp),%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm1,1712(%rsp)
	
# 0 "" 2
# 1318 "../includes/int_double12.0.h" 1
	movapd 1392(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1712(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	shufpd $1,%xmm7,%xmm7
	movapd 1376(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1696(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm7,%xmm0
	movapd %xmm0,1728(%rsp)
	movapd 1392(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1696(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	movapd 1376(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1712(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm0,%xmm7
	movapd %xmm7,1744(%rsp)
# 0 "" 2
# 1003 "../includes/int_double12.0.h" 1
	movapd 1760(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,292(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 292(%rsp)
	je	.L164
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	1760(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L164:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 1760(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 1744(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rbx
	movq	256(%rsp), %rdx
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd 1760(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,292(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 292(%rsp)
	je	.L165
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	1760(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L165:
	flds	.LC19(%rip)
	fstpl	1448(%rsp)
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 1760(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 1728(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 1416(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 1408(%rsp)
	movq	%rbx, 1432(%rsp)
	movq	%rdx, 1424(%rsp)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup 1448(%rsp),%xmm0
	movapd 1424(%rsp),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rdx
	movq	256(%rsp), %rax
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup 1448(%rsp),%xmm0
	movapd 1408(%rsp),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rcx
	movq	%rcx, 1464(%rsp)
	movq	256(%rsp), %rcx
	movq	%rcx, 1456(%rsp)
	movq	%rdx, 1480(%rsp)
	movq	%rax, 1472(%rsp)
#APP
# 1048 "../includes/int_double12.0.h" 1
	movapd 160(%rsp),%xmm0
	addpd 1456(%rsp),%xmm0
	movapd %xmm0,624(%rsp)
	movapd 176(%rsp),%xmm1
	addpd 1472(%rsp),%xmm1
	movapd %xmm1,640(%rsp)
	
# 0 "" 2
#NO_APP
	movq	624(%rsp), %rax
	movq	%rax, 160(%rsp)
	movq	632(%rsp), %rax
	movq	%rax, 168(%rsp)
	movq	640(%rsp), %rax
	movq	%rax, 176(%rsp)
	movq	648(%rsp), %rax
	movq	%rax, 184(%rsp)
	movq	64(%rsp), %rax
	movq	%rax, 128(%rsp)
	movq	72(%rsp), %rax
	movq	%rax, 136(%rsp)
	movq	80(%rsp), %rax
	movq	%rax, 144(%rsp)
	movq	88(%rsp), %rax
	movq	%rax, 152(%rsp)
	movl	$bernoulli+16, %ebx
	movl	$3, %edx
	movl	$d_zero, %r8d
	movl	$d_zero+8, %ebp
	fldz
.L170:
#APP
# 1318 "../includes/int_double12.0.h" 1
	movapd 144(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 80(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	shufpd $1,%xmm7,%xmm7
	movapd 128(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 64(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm7,%xmm0
	movapd %xmm0,1344(%rsp)
	movapd 144(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 64(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	movapd 128(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 80(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm0,%xmm7
	movapd %xmm7,1360(%rsp)
# 0 "" 2
# 1318 "../includes/int_double12.0.h" 1
	movapd 1360(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 80(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	shufpd $1,%xmm7,%xmm7
	movapd 1344(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 64(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm7,%xmm0
	movapd %xmm0,592(%rsp)
	movapd 1360(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 64(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	movapd 1344(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 80(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm0,%xmm7
	movapd %xmm7,608(%rsp)
# 0 "" 2
#NO_APP
	movq	592(%rsp), %rax
	movq	%rax, 128(%rsp)
	movq	600(%rsp), %rax
	movq	%rax, 136(%rsp)
	movq	608(%rsp), %rax
	movq	%rax, 144(%rsp)
	movq	616(%rsp), %rax
	movq	%rax, 152(%rsp)
	movq	(%rbx), %rax
	movq	%rax, 1232(%rsp)
	movq	8(%rbx), %rax
	movq	%rax, 1240(%rsp)
	movq	(%r8), %rax
	movq	%rax, 1248(%rsp)
	movq	0(%rbp), %rax
	movq	%rax, 1256(%rsp)
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd 144(%rsp),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 1784(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 1776(%rsp)
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd 128(%rsp),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 1800(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 1792(%rsp)
#APP
# 392 "../includes/int_double12.0.h" 1
	movapd 1792(%rsp),%XMM0
	addpd 1776(%rsp),%XMM0
	movapd %XMM0,256(%rsp)
	
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 1880(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 1872(%rsp)
#APP
# 1033 "../includes/int_double12.0.h" 1
	movapd 128(%rsp),%xmm0
	movapd %xmm0,1808(%rsp)
	movapd 144(%rsp),%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm1,1824(%rsp)
	
# 0 "" 2
# 1318 "../includes/int_double12.0.h" 1
	movapd 1248(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1824(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	shufpd $1,%xmm7,%xmm7
	movapd 1232(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1808(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm7,%xmm0
	movapd %xmm0,1840(%rsp)
	movapd 1248(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1808(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	movapd 1232(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1824(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm0,%xmm7
	movapd %xmm7,1856(%rsp)
# 0 "" 2
# 1003 "../includes/int_double12.0.h" 1
	movapd 1872(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,292(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 292(%rsp)
	je	.L166
	fstp	%st(0)
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	1872(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L166:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 1872(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 1856(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rdi
	movq	256(%rsp), %rsi
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd 1872(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,292(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 292(%rsp)
	je	.L167
	fstp	%st(0)
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	1872(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L167:
	leal	1(%rdx), %eax
	imull	%edx, %eax
	movq	%rax, 8(%rsp)
	fildll	8(%rsp)
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 1872(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 1840(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rcx
	movq	256(%rsp), %rax
	fld	%st(1)
	fxch	%st(1)
	fucomi	%st(1), %st
	fstp	%st(1)
	jbe	.L192
	fstpl	1304(%rsp)
	movq	%rcx, 1272(%rsp)
	movq	%rax, 1264(%rsp)
	movq	%rdi, 1288(%rsp)
	movq	%rsi, 1280(%rsp)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup 1304(%rsp),%xmm0
	movapd 1280(%rsp),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rcx
	movq	256(%rsp), %rax
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup 1304(%rsp),%xmm0
	movapd 1264(%rsp),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rsi
	movq	%rsi, 1320(%rsp)
	movq	256(%rsp), %rsi
	movq	%rsi, 1312(%rsp)
	movq	%rcx, 1336(%rsp)
	movq	%rax, 1328(%rsp)
#APP
# 1048 "../includes/int_double12.0.h" 1
	movapd 160(%rsp),%xmm0
	addpd 1312(%rsp),%xmm0
	movapd %xmm0,560(%rsp)
	movapd 176(%rsp),%xmm1
	addpd 1328(%rsp),%xmm1
	movapd %xmm1,576(%rsp)
	
# 0 "" 2
#NO_APP
	movq	560(%rsp), %rax
	movq	%rax, 160(%rsp)
	movq	568(%rsp), %rax
	movq	%rax, 168(%rsp)
	movq	576(%rsp), %rax
	movq	%rax, 176(%rsp)
	movq	584(%rsp), %rax
	movq	%rax, 184(%rsp)
	addq	$16, %rbx
	addl	$2, %edx
	cmpl	$15, %edx
	jne	.L170
	fstp	%st(0)
	jmp	.L196
.L192:
	fstp	%st(0)
	fstp	%st(0)
	movl	$.LC49, %edi
	call	puts
	movl	$1, %edi
	call	exit
.L196:
#APP
# 1048 "../includes/int_double12.0.h" 1
	movapd 160(%rsp),%xmm0
	addpd ln_gamma_err(%rip),%xmm0
	movapd %xmm0,192(%rsp)
	movapd 176(%rsp),%xmm1
	addpd ln_gamma_err+16(%rip),%xmm1
	movapd %xmm1,208(%rsp)
	
# 0 "" 2
#NO_APP
	jmp	.L172
.L191:
	movq	c_zero(%rip), %rax
	movq	%rax, 96(%rsp)
	movq	c_zero+8(%rip), %rax
	movq	%rax, 104(%rsp)
	movq	c_zero+16(%rip), %rax
	movq	%rax, 112(%rsp)
	movq	c_zero+24(%rip), %rax
	movq	%rax, 120(%rsp)
	movl	$0, %ebx
	leaq	1200(%rsp), %rbp
	leaq	1168(%rsp), %r12
.L173:
	mov	%ebx, %eax
	movq	%rax, 8(%rsp)
	fildll	8(%rsp)
	fstpl	1160(%rsp)
#APP
# 1091 "../includes/int_double12.0.h" 1
	movddup 1160(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	addpd 64(%rsp),%xmm0
	movapd %xmm0,1168(%rsp)
	movapd 80(%rsp),%xmm1
	movapd %xmm1,1184(%rsp)
	
# 0 "" 2
#NO_APP
	movq	%r12, %rsi
	movq	%rbp, %rdi
	call	_Z3logRK11int_complex
#APP
# 1048 "../includes/int_double12.0.h" 1
	movapd 96(%rsp),%xmm0
	addpd 1200(%rsp),%xmm0
	movapd %xmm0,528(%rsp)
	movapd 112(%rsp),%xmm1
	addpd 1216(%rsp),%xmm1
	movapd %xmm1,544(%rsp)
	
# 0 "" 2
#NO_APP
	movq	528(%rsp), %rax
	movq	%rax, 96(%rsp)
	movq	536(%rsp), %rax
	movq	%rax, 104(%rsp)
	movq	544(%rsp), %rax
	movq	%rax, 112(%rsp)
	movq	552(%rsp), %rax
	movq	%rax, 120(%rsp)
	addl	$1, %ebx
	cmpl	$10, %ebx
	jne	.L173
	flds	.LC52(%rip)
	fstpl	1152(%rsp)
#APP
# 1091 "../includes/int_double12.0.h" 1
	movddup 1152(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	addpd 64(%rsp),%xmm0
	movapd %xmm0,496(%rsp)
	movapd 80(%rsp),%xmm1
	movapd %xmm1,512(%rsp)
	
# 0 "" 2
#NO_APP
	movq	496(%rsp), %rax
	movq	%rax, 64(%rsp)
	movq	504(%rsp), %rax
	movq	%rax, 72(%rsp)
	movq	512(%rsp), %rax
	movq	%rax, 80(%rsp)
	movq	520(%rsp), %rax
	movq	%rax, 88(%rsp)
	flds	.LC3(%rip)
	fstpl	1128(%rsp)
#APP
# 561 "../includes/int_double12.0.h" 1
	movddup 1128(%rsp),%xmm0
	movapd %xmm0,%xmm1
	xorpd d_neg_zero(%rip),%xmm1
	movapd d_ln_two_pi(%rip),%xmm2
	mulpd %xmm2,%xmm0
	mulpd %xmm2,%xmm1
	shufpd $1,%xmm1,%xmm1
	minpd %xmm1,%xmm0
	movapd %xmm0,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 1144(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 1136(%rsp)
	leaq	1024(%rsp), %rdi
	leaq	64(%rsp), %rsi
	call	_Z3logRK11int_complex
	flds	.LC3(%rip)
	fstpl	984(%rsp)
#APP
# 1187 "../includes/int_double12.0.h" 1
	movddup 984(%rsp),%xmm0
	movapd 64(%rsp),%xmm1
	addsubpd %xmm0,%xmm1
	movapd %xmm1,992(%rsp)
	movapd 80(%rsp),%xmm1
	movapd %xmm1,1008(%rsp)
	
# 0 "" 2
# 1318 "../includes/int_double12.0.h" 1
	movapd 1008(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1040(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	shufpd $1,%xmm7,%xmm7
	movapd 992(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1024(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm7,%xmm0
	movapd %xmm0,1056(%rsp)
	movapd 1008(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1024(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	movapd 992(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1040(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm0,%xmm7
	movapd %xmm7,1072(%rsp)
# 0 "" 2
# 1157 "../includes/int_double12.0.h" 1
	movapd 64(%rsp),%xmm0
	shufpd $1,%xmm0,%xmm0
	addpd 1056(%rsp),%xmm0
	movapd %xmm0,1088(%rsp)
	movapd 80(%rsp),%xmm1
	shufpd $1,%xmm1,%xmm1
	addpd 1072(%rsp),%xmm1
	movapd %xmm1,1104(%rsp)
	
# 0 "" 2
# 1062 "../includes/int_double12.0.h" 1
	movapd 1088(%rsp),%xmm0
	addpd 1136(%rsp),%xmm0
	movapd %xmm0,464(%rsp)
	movapd 1104(%rsp),%xmm1
	movapd %xmm1,480(%rsp)
	
# 0 "" 2
#NO_APP
	movq	464(%rsp), %rax
	movq	%rax, 160(%rsp)
	movq	472(%rsp), %rax
	movq	%rax, 168(%rsp)
	movq	480(%rsp), %rax
	movq	%rax, 176(%rsp)
	movq	488(%rsp), %rax
	movq	%rax, 184(%rsp)
	movq	bernoulli(%rip), %rax
	movq	%rax, 864(%rsp)
	movq	bernoulli+8(%rip), %rax
	movq	%rax, 872(%rsp)
	movq	d_zero(%rip), %rax
	movq	%rax, 880(%rsp)
	movq	d_zero+8(%rip), %rax
	movq	%rax, 888(%rsp)
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd 80(%rsp),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 1896(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 1888(%rsp)
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd 64(%rsp),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 1912(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 1904(%rsp)
#APP
# 392 "../includes/int_double12.0.h" 1
	movapd 1904(%rsp),%XMM0
	addpd 1888(%rsp),%XMM0
	movapd %XMM0,256(%rsp)
	
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 1992(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 1984(%rsp)
#APP
# 1033 "../includes/int_double12.0.h" 1
	movapd 64(%rsp),%xmm0
	movapd %xmm0,1920(%rsp)
	movapd 80(%rsp),%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm1,1936(%rsp)
	
# 0 "" 2
# 1318 "../includes/int_double12.0.h" 1
	movapd 880(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1936(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	shufpd $1,%xmm7,%xmm7
	movapd 864(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1920(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm7,%xmm0
	movapd %xmm0,1952(%rsp)
	movapd 880(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1920(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	movapd 864(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1936(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm0,%xmm7
	movapd %xmm7,1968(%rsp)
# 0 "" 2
# 1003 "../includes/int_double12.0.h" 1
	movapd 1984(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,292(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 292(%rsp)
	je	.L174
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	1984(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L174:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 1984(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 1968(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rbx
	movq	256(%rsp), %rdx
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd 1984(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,292(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 292(%rsp)
	je	.L175
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	1984(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L175:
	flds	.LC19(%rip)
	fstpl	936(%rsp)
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 1984(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 1952(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 904(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 896(%rsp)
	movq	%rbx, 920(%rsp)
	movq	%rdx, 912(%rsp)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup 936(%rsp),%xmm0
	movapd 912(%rsp),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rdx
	movq	256(%rsp), %rax
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup 936(%rsp),%xmm0
	movapd 896(%rsp),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rcx
	movq	%rcx, 952(%rsp)
	movq	256(%rsp), %rcx
	movq	%rcx, 944(%rsp)
	movq	%rdx, 968(%rsp)
	movq	%rax, 960(%rsp)
#APP
# 1048 "../includes/int_double12.0.h" 1
	movapd 160(%rsp),%xmm0
	addpd 944(%rsp),%xmm0
	movapd %xmm0,432(%rsp)
	movapd 176(%rsp),%xmm1
	addpd 960(%rsp),%xmm1
	movapd %xmm1,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	432(%rsp), %rax
	movq	%rax, 160(%rsp)
	movq	440(%rsp), %rax
	movq	%rax, 168(%rsp)
	movq	448(%rsp), %rax
	movq	%rax, 176(%rsp)
	movq	456(%rsp), %rax
	movq	%rax, 184(%rsp)
	movq	64(%rsp), %rax
	movq	%rax, 128(%rsp)
	movq	72(%rsp), %rax
	movq	%rax, 136(%rsp)
	movq	80(%rsp), %rax
	movq	%rax, 144(%rsp)
	movq	88(%rsp), %rax
	movq	%rax, 152(%rsp)
	movl	$bernoulli+16, %ebx
	movl	$3, %edx
	movl	$d_zero, %r8d
	movl	$d_zero+8, %ebp
	fldz
.L180:
#APP
# 1318 "../includes/int_double12.0.h" 1
	movapd 144(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 80(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	shufpd $1,%xmm7,%xmm7
	movapd 128(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 64(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm7,%xmm0
	movapd %xmm0,832(%rsp)
	movapd 144(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 64(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	movapd 128(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 80(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm0,%xmm7
	movapd %xmm7,848(%rsp)
# 0 "" 2
# 1318 "../includes/int_double12.0.h" 1
	movapd 848(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 80(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	shufpd $1,%xmm7,%xmm7
	movapd 832(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 64(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm7,%xmm0
	movapd %xmm0,400(%rsp)
	movapd 848(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 64(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	movapd 832(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 80(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm0,%xmm7
	movapd %xmm7,416(%rsp)
# 0 "" 2
#NO_APP
	movq	400(%rsp), %rax
	movq	%rax, 128(%rsp)
	movq	408(%rsp), %rax
	movq	%rax, 136(%rsp)
	movq	416(%rsp), %rax
	movq	%rax, 144(%rsp)
	movq	424(%rsp), %rax
	movq	%rax, 152(%rsp)
	movq	(%rbx), %rax
	movq	%rax, 720(%rsp)
	movq	8(%rbx), %rax
	movq	%rax, 728(%rsp)
	movq	(%r8), %rax
	movq	%rax, 736(%rsp)
	movq	0(%rbp), %rax
	movq	%rax, 744(%rsp)
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd 144(%rsp),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 2008(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 2000(%rsp)
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd 128(%rsp),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 2024(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 2016(%rsp)
#APP
# 392 "../includes/int_double12.0.h" 1
	movapd 2016(%rsp),%XMM0
	addpd 2000(%rsp),%XMM0
	movapd %XMM0,256(%rsp)
	
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 2104(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 2096(%rsp)
#APP
# 1033 "../includes/int_double12.0.h" 1
	movapd 128(%rsp),%xmm0
	movapd %xmm0,2032(%rsp)
	movapd 144(%rsp),%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm1,2048(%rsp)
	
# 0 "" 2
# 1318 "../includes/int_double12.0.h" 1
	movapd 736(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 2048(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	shufpd $1,%xmm7,%xmm7
	movapd 720(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 2032(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm7,%xmm0
	movapd %xmm0,2064(%rsp)
	movapd 736(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 2032(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	movapd 720(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 2048(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm0,%xmm7
	movapd %xmm7,2080(%rsp)
# 0 "" 2
# 1003 "../includes/int_double12.0.h" 1
	movapd 2096(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,292(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 292(%rsp)
	je	.L176
	fstp	%st(0)
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	2096(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L176:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 2096(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 2080(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rdi
	movq	256(%rsp), %rsi
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd 2096(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,292(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 292(%rsp)
	je	.L177
	fstp	%st(0)
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	2096(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L177:
	leal	1(%rdx), %eax
	imull	%edx, %eax
	movq	%rax, 8(%rsp)
	fildll	8(%rsp)
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 2096(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 2064(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rcx
	movq	256(%rsp), %rax
	fld	%st(1)
	fxch	%st(1)
	fucomi	%st(1), %st
	fstp	%st(1)
	jbe	.L194
	fstpl	792(%rsp)
	movq	%rcx, 760(%rsp)
	movq	%rax, 752(%rsp)
	movq	%rdi, 776(%rsp)
	movq	%rsi, 768(%rsp)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup 792(%rsp),%xmm0
	movapd 768(%rsp),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rcx
	movq	256(%rsp), %rax
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup 792(%rsp),%xmm0
	movapd 752(%rsp),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rsi
	movq	%rsi, 808(%rsp)
	movq	256(%rsp), %rsi
	movq	%rsi, 800(%rsp)
	movq	%rcx, 824(%rsp)
	movq	%rax, 816(%rsp)
#APP
# 1048 "../includes/int_double12.0.h" 1
	movapd 160(%rsp),%xmm0
	addpd 800(%rsp),%xmm0
	movapd %xmm0,368(%rsp)
	movapd 176(%rsp),%xmm1
	addpd 816(%rsp),%xmm1
	movapd %xmm1,384(%rsp)
	
# 0 "" 2
#NO_APP
	movq	368(%rsp), %rax
	movq	%rax, 160(%rsp)
	movq	376(%rsp), %rax
	movq	%rax, 168(%rsp)
	movq	384(%rsp), %rax
	movq	%rax, 176(%rsp)
	movq	392(%rsp), %rax
	movq	%rax, 184(%rsp)
	addq	$16, %rbx
	addl	$2, %edx
	cmpl	$15, %edx
	jne	.L180
	fstp	%st(0)
	jmp	.L197
.L194:
	fstp	%st(0)
	fstp	%st(0)
	movl	$.LC49, %edi
	call	puts
	movl	$1, %edi
	call	exit
.L197:
#APP
# 1157 "../includes/int_double12.0.h" 1
	movapd 96(%rsp),%xmm0
	shufpd $1,%xmm0,%xmm0
	addpd 160(%rsp),%xmm0
	movapd %xmm0,688(%rsp)
	movapd 112(%rsp),%xmm1
	shufpd $1,%xmm1,%xmm1
	addpd 176(%rsp),%xmm1
	movapd %xmm1,704(%rsp)
	
# 0 "" 2
# 1048 "../includes/int_double12.0.h" 1
	movapd 688(%rsp),%xmm0
	addpd ln_gamma_err_2(%rip),%xmm0
	movapd %xmm0,192(%rsp)
	movapd 704(%rsp),%xmm1
	addpd ln_gamma_err_2+16(%rip),%xmm1
	movapd %xmm1,208(%rsp)
	
# 0 "" 2
#NO_APP
.L172:
#APP
# 524 "../includes/int_double12.0.h" 1
	movapd 240(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd d_ln_pi(%rip),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,256(%rsp)
	
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 328(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 320(%rsp)
#APP
# 454 "../includes/int_double12.0.h" 1
	movapd 320(%rsp),%xmm0
	shufpd $1,%xmm0,%xmm0
	addpd 208(%rsp),%xmm0
	movapd %xmm0,256(%rsp)
	
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rbx
	movq	256(%rsp), %rdx
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd d_pi(%rip),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,292(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 292(%rsp)
	je	.L182
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	movl	$d_pi, %edi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L182:
	movq	%rbx, 312(%rsp)
	movq	%rdx, 304(%rsp)
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd d_pi(%rip),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 304(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,256(%rsp)
# 0 "" 2
#NO_APP
	movq	264(%rsp), %rax
	movq	%rax, 280(%rsp)
	movq	256(%rsp), %rax
	movq	%rax, 272(%rsp)
	movq	%r13, 8(%rsp)
	fildll	8(%rsp)
	testq	%r13, %r13
	jns	.L183
	fadds	.LC38(%rip)
.L183:
	fstpl	40(%rsp)
	fldl	40(%rsp)
	fstl	352(%rsp)
	fchs
	fstpl	360(%rsp)
#APP
# 454 "../includes/int_double12.0.h" 1
	movapd 272(%rsp),%xmm0
	shufpd $1,%xmm0,%xmm0
	addpd 352(%rsp),%xmm0
	movapd %xmm0,256(%rsp)
	
# 0 "" 2
#NO_APP
	movsd	256(%rsp), %xmm0
	movsd	264(%rsp), %xmm1
	addq	$2120, %rsp
	.cfi_def_cfa_offset 40
	popq	%rbx
	.cfi_def_cfa_offset 32
	popq	%rbp
	.cfi_def_cfa_offset 24
	popq	%r12
	.cfi_def_cfa_offset 16
	popq	%r13
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE1187:
	.size	_Z1S10int_doublem, .-_Z1S10int_doublem
	.section	.rodata.str1.1
.LC53:
	.string	"Usage:- %s <zeros file>.\n"
.LC56:
	.string	"rb"
	.section	.rodata.str1.8
	.align 8
.LC57:
	.string	"Error opening file %s for binary input. Exiting.\n"
	.align 8
.LC58:
	.string	"Error reading number of iterations from zeros file. Exiting."
	.section	.rodata.str1.1
.LC60:
	.string	"st[0] was 0.0. Exiting."
.LC61:
	.string	"get_zero returned 0. Exiting."
.LC64:
	.string	"0: First Zero %lu is at "
.LC65:
	.string	"z=%lu last_t="
.LC66:
	.string	" t="
	.section	.rodata.str1.8
	.align 8
.LC67:
	.string	"del_t=%lu %lu %lu =[ %20.18e , %20.18e ]\n"
	.section	.rodata.str1.1
.LC68:
	.string	"sum_zeros=%lu %lu %lu\n"
.LC69:
	.string	"1: Last Zero %lu is at "
	.section	.rodata.str1.8
	.align 8
.LC70:
	.string	"2: There is no absolute gap larger than            %20.18e %lu\n"
	.align 8
.LC71:
	.string	"3: There is an absolute as large as                %20.18e %lu\n"
	.align 8
.LC72:
	.string	"4: There is no absolute gap smaller than           %20.18e %lu\n"
	.align 8
.LC73:
	.string	"5: There is an absolute gap as small as            %20.18e %lu\n"
	.align 8
.LC74:
	.string	"6: There is no relative gap larger than            %20.18e %lu\n"
	.align 8
.LC75:
	.string	"7: There is a relative as large as                 %20.18e %lu\n"
	.align 8
.LC76:
	.string	"8: There is no relative gap smaller than           %20.18e %lu\n"
	.align 8
.LC77:
	.string	"9: There is a relative gap as small as             %20.18e %lu\n"
	.align 8
.LC78:
	.string	"10: There is an S(t) larger than                   %20.18e %lu\n"
	.align 8
.LC79:
	.string	"11: There is no S(t) larger than                   %20.18e %lu\n"
	.align 8
.LC80:
	.string	"12: There is no S(t) smaller than                  %20.18e %lu\n"
	.align 8
.LC81:
	.string	"13: There is an S(t) smaller than                  %20.18e %lu\n"
	.section	.rodata.str1.1
.LC82:
	.string	"Sum 1/gamma         = "
.LC83:
	.string	"Sum 1/rho+1/(1-rho) = "
	.text
.globl main
	.type	main, @function
main:
.LFB1188:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$824, %rsp
	.cfi_def_cfa_offset 880
	movq	%rsi, %rbx
	cmpl	$2, %edi
	je	.L199
	movq	(%rsi), %rsi
	movl	$.LC53, %edi
	movl	$0, %eax
	call	printf
	movl	$0, %edi
	call	exit
.L199:
	call	_Z9_fpu_rnddv
	flds	.LC54(%rip)
	fstl	512(%rsp)
	fstpl	520(%rsp)
	flds	.LC55(%rip)
	fstl	496(%rsp)
	fstpl	504(%rsp)
	addq	$8, %rbx
	movl	$.LC56, %esi
	movq	(%rbx), %rdi
	call	fopen
	movq	%rax, %r13
	testq	%rax, %rax
	jne	.L200
	movq	(%rbx), %rsi
	movl	$.LC57, %edi
	movl	$0, %eax
	call	printf
	movl	$0, %edi
	call	exit
.L200:
	leaq	528(%rsp), %rdi
	movq	%rax, %rcx
	movl	$1, %edx
	movl	$8, %esi
	call	fread
	cmpq	$1, %rax
	je	.L201
	movl	$.LC58, %edi
	call	puts
	movl	$0, %edi
	call	exit
.L201:
	fldz
	fstl	416(%rsp)
	fldz
	fchs
	fstl	424(%rsp)
	fxch	%st(1)
	fstl	400(%rsp)
	fxch	%st(1)
	fstpl	408(%rsp)
	cmpq	$0, 528(%rsp)
	jne	.L202
	fld	%st(0)
	fstl	40(%rsp)
	fstpl	48(%rsp)
	fldl	.LC59(%rip)
	fstl	32(%rsp)
	fstl	56(%rsp)
	fxch	%st(1)
	fstl	88(%rsp)
	fstl	96(%rsp)
	fxch	%st(1)
	fstl	104(%rsp)
	fstl	112(%rsp)
	fxch	%st(1)
	fstl	64(%rsp)
	fstpl	72(%rsp)
	fstl	80(%rsp)
	fstpl	120(%rsp)
	jmp	.L203
.L202:
	fstp	%st(0)
	fldz
	fstl	40(%rsp)
	fstl	48(%rsp)
	fldl	.LC59(%rip)
	fstl	32(%rsp)
	fstl	56(%rsp)
	fxch	%st(1)
	fstl	88(%rsp)
	fstl	96(%rsp)
	fxch	%st(1)
	fstl	104(%rsp)
	fstl	112(%rsp)
	fxch	%st(1)
	fstl	64(%rsp)
	fstpl	72(%rsp)
	fstl	80(%rsp)
	fstpl	120(%rsp)
	movl	$0, %ebx
	movq	$0, 216(%rsp)
	leaq	464(%rsp), %rax
	addq	$8, %rax
	movq	%rax, 224(%rsp)
	leaq	272(%rsp), %r14
.L243:
	movq	$0, 320(%rsp)
	movl	$0, 312(%rsp)
	movq	$0, 304(%rsp)
	movq	%r13, %rcx
	movl	$2, %edx
	movl	$8, %esi
	leaq	480(%rsp), %rdi
	call	fread
	movq	%r13, %rcx
	movl	$1, %edx
	movl	$8, %esi
	leaq	464(%rsp), %rdi
	call	fread
	fldz
	fldl	480(%rsp)
	fucomip	%st(1), %st
	fstp	%st(0)
	jne	.L262
	jp	.L262
	movl	$.LC60, %edi
	call	puts
	movl	$0, %edi
	call	exit
.L262:
	movq	%r13, %rcx
	movl	$1, %edx
	movl	$8, %esi
	movq	224(%rsp), %rdi
	call	fread
	movq	464(%rsp), %rbp
	addq	$1, %rbp
	cmpq	472(%rsp), %rbp
	ja	.L206
.L259:
	movq	%r13, %rsi
	movq	%r14, %rdi
	call	_Z8get_zeroP8_IO_FILE
	cmpq	$0, 272(%rsp)
	jne	.L207
	cmpl	$0, 280(%rsp)
	jne	.L207
	cmpq	$0, 288(%rsp)
	jne	.L207
	movl	$.LC61, %edi
	call	puts
	movl	$0, %edi
	call	exit
.L207:
#APP
# 66 "debug.cpp" 1
	movq 304(%rsp),%r12
	addq 272(%rsp),%r12
	movq %r12,304(%rsp)
	movl 312(%rsp),%r12d
	adcl 280(%rsp),%r12d
	movl %r12d,312(%rsp)
	movq 320(%rsp),%r12
	adcq 288(%rsp),%r12
	movq %r12,320(%rsp)
	
# 0 "" 2
#NO_APP
	fildll	320(%rsp)
	cmpq	$0, 320(%rsp)
	jns	.L208
	fadds	.LC38(%rip)
.L208:
	fstpl	264(%rsp)
	fldl	264(%rsp)
	fstl	368(%rsp)
	fchs
	fstpl	376(%rsp)
	flds	.LC62(%rip)
	fstpl	768(%rsp)
#APP
# 561 "../includes/int_double12.0.h" 1
	movddup 768(%rsp),%xmm0
	movapd %xmm0,%xmm1
	xorpd d_neg_zero(%rip),%xmm1
	movapd 368(%rsp),%xmm2
	mulpd %xmm2,%xmm0
	mulpd %xmm2,%xmm1
	shufpd $1,%xmm1,%xmm1
	minpd %xmm1,%xmm0
	movapd %xmm0,352(%rsp)
# 0 "" 2
#NO_APP
	movq	360(%rsp), %rax
	movq	%rax, 376(%rsp)
	movq	352(%rsp), %rax
	movq	%rax, 368(%rsp)
	mov	312(%rsp), %eax
	movq	%rax, 232(%rsp)
	fildll	232(%rsp)
	fstpl	760(%rsp)
#APP
# 441 "../includes/int_double12.0.h" 1
	movddup 760(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd 368(%rsp),%XMM0
	movapd %XMM0,368(%rsp)
	
# 0 "" 2
#NO_APP
	flds	.LC38(%rip)
	fstl	752(%rsp)
#APP
# 561 "../includes/int_double12.0.h" 1
	movddup 752(%rsp),%xmm0
	movapd %xmm0,%xmm1
	xorpd d_neg_zero(%rip),%xmm1
	movapd 368(%rsp),%xmm2
	mulpd %xmm2,%xmm0
	mulpd %xmm2,%xmm1
	shufpd $1,%xmm1,%xmm1
	minpd %xmm1,%xmm0
	movapd %xmm0,352(%rsp)
# 0 "" 2
#NO_APP
	movq	360(%rsp), %rax
	movq	%rax, 376(%rsp)
	movq	352(%rsp), %rax
	movq	%rax, 368(%rsp)
	fildll	304(%rsp)
	cmpq	$0, 304(%rsp)
	jns	.L264
	faddp	%st, %st(1)
	jmp	.L209
.L264:
	fstp	%st(1)
.L209:
	fstpl	744(%rsp)
#APP
# 441 "../includes/int_double12.0.h" 1
	movddup 744(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd 368(%rsp),%XMM0
	movapd %XMM0,368(%rsp)
	
# 0 "" 2
#NO_APP
	flds	.LC63(%rip)
	fstpl	736(%rsp)
#APP
# 561 "../includes/int_double12.0.h" 1
	movddup 736(%rsp),%xmm0
	movapd %xmm0,%xmm1
	xorpd d_neg_zero(%rip),%xmm1
	movapd 368(%rsp),%xmm2
	mulpd %xmm2,%xmm0
	mulpd %xmm2,%xmm1
	shufpd $1,%xmm1,%xmm1
	minpd %xmm1,%xmm0
	movapd %xmm0,352(%rsp)
# 0 "" 2
#NO_APP
	movq	360(%rsp), %rdx
	movq	352(%rsp), %rax
	movq	%rdx, 376(%rsp)
	movq	%rax, 368(%rsp)
	movq	%rdx, 552(%rsp)
	movq	%rax, 544(%rsp)
#APP
# 392 "../includes/int_double12.0.h" 1
	movapd 544(%rsp),%XMM0
	addpd 512(%rsp),%XMM0
	movapd %XMM0,352(%rsp)
	
# 0 "" 2
#NO_APP
	movq	360(%rsp), %rax
	movq	%rax, 456(%rsp)
	movq	352(%rsp), %rax
	movq	%rax, 448(%rsp)
#APP
# 441 "../includes/int_double12.0.h" 1
	movddup 480(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd 448(%rsp),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	fld1
	fstpl	576(%rsp)
	fld1
	fchs
	fstpl	584(%rsp)
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd 448(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,540(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 540(%rsp)
	je	.L210
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	448(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L210:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 448(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 576(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,336(%rsp)
# 0 "" 2
#NO_APP
	movq	344(%rsp), %rax
	movq	%rax, 568(%rsp)
	movq	336(%rsp), %rax
	movq	%rax, 560(%rsp)
#APP
# 429 "../includes/int_double12.0.h" 1
	movapd 400(%rsp),%XMM0
	addpd 560(%rsp),%XMM0
	movapd %XMM0,400(%rsp)
	
# 0 "" 2
# 961 "../includes/int_double12.0.h" 1
	movapd 448(%rsp),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,352(%rsp)
# 0 "" 2
#NO_APP
	movq	360(%rsp), %rax
	movq	%rax, 632(%rsp)
	movq	352(%rsp), %rax
	movq	%rax, 624(%rsp)
	flds	.LC51(%rip)
	fstpl	648(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 648(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd 624(%rsp),%XMM0
	movapd %XMM0,352(%rsp)
	
# 0 "" 2
#NO_APP
	movq	360(%rsp), %rax
	movq	%rax, 616(%rsp)
	movq	352(%rsp), %rax
	movq	%rax, 608(%rsp)
	fld1
	fstpl	656(%rsp)
	fld1
	fchs
	fstpl	664(%rsp)
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd 608(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,540(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 540(%rsp)
	je	.L211
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	608(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L211:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 608(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 656(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,336(%rsp)
# 0 "" 2
#NO_APP
	movq	344(%rsp), %rax
	movq	%rax, 600(%rsp)
	movq	336(%rsp), %rax
	movq	%rax, 592(%rsp)
#APP
# 429 "../includes/int_double12.0.h" 1
	movapd 416(%rsp),%XMM0
	addpd 592(%rsp),%XMM0
	movapd %XMM0,416(%rsp)
	
# 0 "" 2
#NO_APP
	leaq	-1(%rbp), %rdi
	movsd	448(%rsp), %xmm0
	movsd	456(%rsp), %xmm1
	call	_Z1S10int_doublem
	movsd	%xmm0, 240(%rsp)
	movsd	%xmm1, 248(%rsp)
	movq	240(%rsp), %rax
	movq	%rax, 384(%rsp)
	movq	248(%rsp), %rax
	movq	%rax, 392(%rsp)
	fldl	384(%rsp)
	fldl	48(%rsp)
	fxch	%st(1)
	fucomi	%st(1), %st
	fstp	%st(1)
	jbe	.L265
	movq	%rbp, %r15
	fstpl	48(%rsp)
	jmp	.L212
.L265:
	fstp	%st(0)
.L212:
	fldl	392(%rsp)
	fchs
	fldl	40(%rsp)
	fxch	%st(1)
	fucomi	%st(1), %st
	fstp	%st(1)
	jbe	.L266
	movq	%rbp, 128(%rsp)
	fstpl	40(%rsp)
	jmp	.L214
.L266:
	fstp	%st(0)
.L214:
	fld1
	fstpl	680(%rsp)
#APP
# 467 "../includes/int_double12.0.h" 1
	movddup 680(%rsp),%XMM0
	movapd 384(%rsp),%XMM1
	addsubpd %XMM0,%XMM1
	movapd %XMM1,352(%rsp)
# 0 "" 2
#NO_APP
	fldl	360(%rsp)
	fldl	352(%rsp)
	fxch	%st(1)
	fstl	392(%rsp)
	fxch	%st(1)
	fstl	384(%rsp)
	fldl	56(%rsp)
	fucomip	%st(1), %st
	jbe	.L267
	movq	%rbp, 136(%rsp)
	fstpl	56(%rsp)
	jmp	.L216
.L267:
	fstp	%st(0)
.L216:
	fchs
	fldl	32(%rsp)
	fucomip	%st(1), %st
	jbe	.L268
	movq	%rbp, 144(%rsp)
	fstpl	32(%rsp)
	jmp	.L218
.L268:
	fstp	%st(0)
.L218:
	testl	%ebx, %ebx
	jne	.L220
	movq	%rbp, %rsi
	movl	$.LC64, %edi
	movl	$0, %eax
	call	printf
	movl	$32, %edi
	call	putchar
	leaq	448(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movb	$1, %bl
	jmp	.L221
.L220:
	fildll	288(%rsp)
	cmpq	$0, 288(%rsp)
	jns	.L222
	fadds	.LC38(%rip)
.L222:
	fstpl	264(%rsp)
	fldl	264(%rsp)
	fstl	352(%rsp)
	fchs
	fstpl	360(%rsp)
	flds	.LC62(%rip)
	fstpl	808(%rsp)
#APP
# 561 "../includes/int_double12.0.h" 1
	movddup 808(%rsp),%xmm0
	movapd %xmm0,%xmm1
	xorpd d_neg_zero(%rip),%xmm1
	movapd 352(%rsp),%xmm2
	mulpd %xmm2,%xmm0
	mulpd %xmm2,%xmm1
	shufpd $1,%xmm1,%xmm1
	minpd %xmm1,%xmm0
	movapd %xmm0,336(%rsp)
# 0 "" 2
#NO_APP
	movq	344(%rsp), %rax
	movq	%rax, 360(%rsp)
	movq	336(%rsp), %rax
	movq	%rax, 352(%rsp)
	mov	280(%rsp), %eax
	movq	%rax, 232(%rsp)
	fildll	232(%rsp)
	fstpl	800(%rsp)
#APP
# 441 "../includes/int_double12.0.h" 1
	movddup 800(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd 352(%rsp),%XMM0
	movapd %XMM0,352(%rsp)
	
# 0 "" 2
#NO_APP
	flds	.LC38(%rip)
	fstl	792(%rsp)
#APP
# 561 "../includes/int_double12.0.h" 1
	movddup 792(%rsp),%xmm0
	movapd %xmm0,%xmm1
	xorpd d_neg_zero(%rip),%xmm1
	movapd 352(%rsp),%xmm2
	mulpd %xmm2,%xmm0
	mulpd %xmm2,%xmm1
	shufpd $1,%xmm1,%xmm1
	minpd %xmm1,%xmm0
	movapd %xmm0,336(%rsp)
# 0 "" 2
#NO_APP
	movq	344(%rsp), %rax
	movq	%rax, 360(%rsp)
	movq	336(%rsp), %rax
	movq	%rax, 352(%rsp)
	fildll	272(%rsp)
	cmpq	$0, 272(%rsp)
	jns	.L269
	faddp	%st, %st(1)
	jmp	.L223
.L269:
	fstp	%st(1)
.L223:
	fstpl	784(%rsp)
#APP
# 441 "../includes/int_double12.0.h" 1
	movddup 784(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd 352(%rsp),%XMM0
	movapd %XMM0,352(%rsp)
	
# 0 "" 2
#NO_APP
	flds	.LC63(%rip)
	fstpl	776(%rsp)
#APP
# 561 "../includes/int_double12.0.h" 1
	movddup 776(%rsp),%xmm0
	movapd %xmm0,%xmm1
	xorpd d_neg_zero(%rip),%xmm1
	movapd 352(%rsp),%xmm2
	mulpd %xmm2,%xmm0
	mulpd %xmm2,%xmm1
	shufpd $1,%xmm1,%xmm1
	minpd %xmm1,%xmm0
	movapd %xmm0,336(%rsp)
# 0 "" 2
#NO_APP
	movq	344(%rsp), %rdx
	movq	336(%rsp), %rax
	movq	%rdx, 360(%rsp)
	movq	%rax, 352(%rsp)
	movq	%rdx, 696(%rsp)
	movq	%rax, 688(%rsp)
#APP
# 392 "../includes/int_double12.0.h" 1
	movapd 688(%rsp),%XMM0
	addpd 496(%rsp),%XMM0
	movapd %XMM0,336(%rsp)
	
# 0 "" 2
#NO_APP
	fldl	344(%rsp)
	fldl	336(%rsp)
	fxch	%st(1)
	fstl	376(%rsp)
	fxch	%st(1)
	fstl	368(%rsp)
	fldl	72(%rsp)
	fxch	%st(1)
	fucomi	%st(1), %st
	fstp	%st(1)
	jbe	.L270
	movq	%rbp, 152(%rsp)
	fstl	72(%rsp)
	fxch	%st(1)
	jmp	.L224
.L270:
	fxch	%st(1)
.L224:
	fchs
	fldl	64(%rsp)
	fxch	%st(1)
	fucomi	%st(1), %st
	fstp	%st(1)
	jbe	.L271
	movq	%rbp, 160(%rsp)
	fstpl	64(%rsp)
	jmp	.L226
.L271:
	fstp	%st(0)
.L226:
	fldl	120(%rsp)
	fucomip	%st(1), %st
	fstp	%st(0)
	jbe	.L228
	movq	%rbp, %rsi
	movl	$.LC65, %edi
	movl	$0, %eax
	call	printf
	leaq	432(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$.LC66, %edi
	movl	$0, %eax
	call	printf
	movl	$32, %edi
	call	putchar
	leaq	448(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	fldl	376(%rsp)
	fchs
	fstpl	8(%rsp)
	movsd	8(%rsp), %xmm1
	movsd	368(%rsp), %xmm0
	movq	288(%rsp), %rcx
	movl	280(%rsp), %edx
	movq	272(%rsp), %rsi
	movl	$.LC67, %edi
	movl	$2, %eax
	call	printf
	movq	320(%rsp), %rcx
	movl	312(%rsp), %edx
	movq	304(%rsp), %rsi
	movl	$.LC68, %edi
	movl	$0, %eax
	call	printf
	fldl	368(%rsp)
	fstpl	120(%rsp)
	movq	%rbp, 208(%rsp)
.L228:
	fldl	376(%rsp)
	fchs
	fldl	80(%rsp)
	fucomip	%st(1), %st
	jbe	.L272
	movq	%rbp, 168(%rsp)
	fstpl	80(%rsp)
	jmp	.L230
.L272:
	fstp	%st(0)
.L230:
	fldl	440(%rsp)
	fchs
	fstpl	24(%rsp)
	movsd	24(%rsp), %xmm0
	call	log_ru
	movsd	%xmm0, 24(%rsp)
	movsd	432(%rsp), %xmm0
	call	log_rd
	movsd	%xmm0, 8(%rsp)
	fldl	8(%rsp)
	fldl	24(%rsp)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	jbe	.L263
	movapd	%xmm0, %xmm1
	movsd	24(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L263:
	fldl	24(%rsp)
	fchs
	fstpl	728(%rsp)
	movsd	%xmm0, 720(%rsp)
#APP
# 524 "../includes/int_double12.0.h" 1
	movapd 720(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 368(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,336(%rsp)
	
# 0 "" 2
#NO_APP
	movq	344(%rsp), %rax
	movq	%rax, 712(%rsp)
	movq	336(%rsp), %rax
	movq	%rax, 704(%rsp)
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd d_two_pi(%rip),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,540(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 540(%rsp)
	je	.L234
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	movl	$d_two_pi, %edi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L234:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd d_two_pi(%rip),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 704(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,336(%rsp)
# 0 "" 2
#NO_APP
	fldl	344(%rsp)
	fldl	336(%rsp)
	fldl	96(%rsp)
	fxch	%st(1)
	fucomi	%st(1), %st
	fstp	%st(1)
	jbe	.L273
	movq	%rbp, 176(%rsp)
	fstl	96(%rsp)
	fxch	%st(1)
	jmp	.L235
.L273:
	fxch	%st(1)
.L235:
	fchs
	fldl	88(%rsp)
	fxch	%st(1)
	fucomi	%st(1), %st
	fstp	%st(1)
	jbe	.L237
	movq	%rbp, 184(%rsp)
	fstl	88(%rsp)
.L237:
	fldl	112(%rsp)
	fucomip	%st(2), %st
	jbe	.L274
	fxch	%st(1)
	movq	%rbp, 192(%rsp)
	fstpl	112(%rsp)
	jmp	.L239
.L274:
	fstp	%st(1)
.L239:
	fldl	104(%rsp)
	fucomip	%st(1), %st
	jbe	.L275
	movq	%rbp, 200(%rsp)
	fstpl	104(%rsp)
	jmp	.L221
.L275:
	fstp	%st(0)
.L221:
	movq	448(%rsp), %rax
	movq	%rax, 432(%rsp)
	movq	456(%rsp), %rax
	movq	%rax, 440(%rsp)
	addq	$1, %rbp
	cmpq	%rbp, 472(%rsp)
	jae	.L259
.L206:
	addq	$1, 216(%rsp)
	movq	216(%rsp), %rax
	cmpq	%rax, 528(%rsp)
	ja	.L243
.L203:
	leaq	-1(%rbp), %rsi
	movl	$.LC69, %edi
	movl	$0, %eax
	call	printf
	movl	$32, %edi
	call	putchar
	leaq	448(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movq	160(%rsp), %rsi
	movsd	64(%rsp), %xmm0
	movl	$.LC70, %edi
	movl	$1, %eax
	call	printf
	movq	152(%rsp), %rsi
	movsd	72(%rsp), %xmm0
	movl	$.LC71, %edi
	movl	$1, %eax
	call	printf
	movq	208(%rsp), %rsi
	movsd	120(%rsp), %xmm0
	movl	$.LC72, %edi
	movl	$1, %eax
	call	printf
	movq	168(%rsp), %rsi
	movsd	80(%rsp), %xmm0
	movl	$.LC73, %edi
	movl	$1, %eax
	call	printf
	movq	184(%rsp), %rsi
	movsd	88(%rsp), %xmm0
	movl	$.LC74, %edi
	movl	$1, %eax
	call	printf
	movq	176(%rsp), %rsi
	movsd	96(%rsp), %xmm0
	movl	$.LC75, %edi
	movl	$1, %eax
	call	printf
	movq	192(%rsp), %rsi
	movsd	112(%rsp), %xmm0
	movl	$.LC76, %edi
	movl	$1, %eax
	call	printf
	movq	200(%rsp), %rsi
	movsd	104(%rsp), %xmm0
	movl	$.LC77, %edi
	movl	$1, %eax
	call	printf
	movq	%r15, %rsi
	movsd	48(%rsp), %xmm0
	movl	$.LC78, %edi
	movl	$1, %eax
	call	printf
	movq	128(%rsp), %rsi
	movsd	40(%rsp), %xmm0
	movl	$.LC79, %edi
	movl	$1, %eax
	call	printf
	movq	136(%rsp), %rsi
	movsd	56(%rsp), %xmm0
	movl	$.LC80, %edi
	movl	$1, %eax
	call	printf
	movq	144(%rsp), %rsi
	movsd	32(%rsp), %xmm0
	movl	$.LC81, %edi
	movl	$1, %eax
	call	printf
	movl	$.LC82, %edi
	movl	$0, %eax
	call	printf
	movl	$32, %edi
	call	putchar
	leaq	400(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$.LC83, %edi
	movl	$0, %eax
	call	printf
	movl	$32, %edi
	call	putchar
	leaq	416(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movq	%r13, %rdi
	call	fclose
	movl	$0, %eax
	addq	$824, %rsp
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE1188:
	.size	main, .-main
	.section	.rodata.str1.8
	.align 8
.LC84:
	.string	"s must be exact and in first quadrant"
	.section	.rodata.str1.1
.LC85:
	.string	"Exiting."
	.text
.globl _Z7hurwitzRK11int_complexRK10int_double
	.type	_Z7hurwitzRK11int_complexRK10int_double, @function
_Z7hurwitzRK11int_complexRK10int_double:
.LFB1163:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$1560, %rsp
	.cfi_def_cfa_offset 1616
	movq	%rdi, %rbx
	movq	%rsi, %rbp
	movq	%rdx, %r13
	fldl	(%rsi)
	fld	%st(0)
	faddl	8(%rsi)
	fldz
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	jne	.L385
	jp	.L386
	fldl	16(%rsi)
	fld	%st(0)
	faddl	24(%rsi)
	fldz
	fxch	%st(1)
	fucomip	%st(1), %st
	jne	.L387
	jp	.L388
	fucomi	%st(2), %st
	fstp	%st(2)
	jae	.L389
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	.p2align 4,,3
	jbe	.L372
	.p2align 4,,3
	jmp	.L277
.L385:
	fstp	%st(0)
	.p2align 4,,5
	jmp	.L277
.L386:
	fstp	%st(0)
	.p2align 4,,7
	jmp	.L277
.L387:
	fstp	%st(0)
	fstp	%st(0)
	fstp	%st(0)
	jmp	.L277
.L388:
	fstp	%st(0)
	fstp	%st(0)
	fstp	%st(0)
	jmp	.L277
.L389:
	fstp	%st(0)
	fstp	%st(0)
.L277:
	movl	$.LC84, %edi
	movl	$0, %eax
	call	printf
	movl	$32, %edi
	call	putchar
	movq	%rbp, %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$43, %edi
	call	putchar
	leaq	16(%rbp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$105, %edi
	call	putchar
	movl	$10, %edi
	call	putchar
	movl	$.LC85, %edi
	call	puts
	movl	$1, %edi
	call	exit
.L372:
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd 16(%rsi),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,464(%rsp)
# 0 "" 2
#NO_APP
	movq	472(%rsp), %rax
	movq	%rax, 1176(%rsp)
	movq	464(%rsp), %rax
	movq	%rax, 1168(%rsp)
	movq	%rsi, 48(%rsp)
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd (%rsi),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,464(%rsp)
# 0 "" 2
#NO_APP
	movq	472(%rsp), %rax
	movq	%rax, 1160(%rsp)
	movq	464(%rsp), %rax
	movq	%rax, 1152(%rsp)
#APP
# 392 "../includes/int_double12.0.h" 1
	movapd 1152(%rsp),%XMM0
	addpd 1168(%rsp),%XMM0
	movapd %XMM0,464(%rsp)
	
# 0 "" 2
#NO_APP
	movq	464(%rsp), %rax
	movq	%rax, 536(%rsp)
	movq	472(%rsp), %rax
	movq	%rax, 544(%rsp)
#APP
# 979 "../includes/int_double12.0.h" 1
	fldl 536(%rsp)
	fsqrt
	fstpl 552(%rsp)
	fldl 544(%rsp)
	fchs
	fsqrt
	fstpl 560(%rsp)
	
# 0 "" 2
#NO_APP
	fldl	560(%rsp)
	fldl	552(%rsp)
	fucomi	%st(1), %st
	jbe	.L373
	fstpl	8(%rsp)
	movsd	8(%rsp), %xmm1
	fstpl	8(%rsp)
	movsd	8(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L373:
	fxch	%st(1)
	fchs
	fstpl	472(%rsp)
	fstpl	464(%rsp)
#APP
# 392 "../includes/int_double12.0.h" 1
	movapd 464(%rsp),%XMM0
	addpd delta_blow(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 488(%rsp)
	movq	448(%rsp), %rax
	movq	%rax, 480(%rsp)
	movq	(%rsi), %rax
	movq	%rax, 416(%rsp)
	movq	8(%rsi), %rax
	movq	%rax, 424(%rsp)
	movq	16(%rsi), %rax
	movq	%rax, 432(%rsp)
	movq	24(%rsi), %rax
	movq	%rax, 440(%rsp)
	movq	c_zero(%rip), %rax
	movq	%rax, (%rdi)
	movq	c_zero+8(%rip), %rax
	movq	%rax, 8(%rdi)
	movq	c_zero+16(%rip), %rax
	movq	%rax, 16(%rdi)
	movq	c_zero+24(%rip), %rax
	movq	%rax, 24(%rdi)
	cmpb	$0, h_bernoulli_initialised(%rip)
	jne	.L282
	call	_Z15set_h_bernoulliv
	movb	$1, h_bernoulli_initialised(%rip)
.L282:
	fldl	488(%rsp)
	fchs
	flds	.LC86(%rip)
	fxch	%st(1)
	fucomi	%st(1), %st
	fstp	%st(1)
	jbe	.L374
	fstpl	8(%rsp)
	movsd	8(%rsp), %xmm0
	call	ceil
	cvttsd2siq	%xmm0, %r14
	mov	%r14d, %eax
	movq	%rax, 56(%rsp)
	fildll	56(%rsp)
	fstpl	568(%rsp)
#APP
# 405 "../includes/int_double12.0.h" 1
	movddup 568(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd 0(%r13),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 504(%rsp)
	movq	448(%rsp), %rax
	movq	%rax, 496(%rsp)
	fld1
	fstpl	616(%rsp)
#APP
# 1140 "../includes/int_double12.0.h" 1
	movapd 0(%rbp),%xmm0
	shufpd $1,%xmm0,%xmm0
	movapd %xmm0,624(%rsp)
	movapd 16(%rbp),%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm1,640(%rsp)
	
# 0 "" 2
# 1091 "../includes/int_double12.0.h" 1
	movddup 616(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	addpd 624(%rsp),%xmm0
	movapd %xmm0,576(%rsp)
	movapd 640(%rsp),%xmm1
	movapd %xmm1,592(%rsp)
	
# 0 "" 2
#NO_APP
	leaq	64(%rsp), %rdi
	leaq	576(%rsp), %rdx
	leaq	496(%rsp), %rsi
	call	_Z3powRK10int_doubleRK11int_complex
	movq	64(%rsp), %rax
	movq	%rax, 384(%rsp)
	movq	72(%rsp), %rax
	movq	%rax, 392(%rsp)
	movq	80(%rsp), %rax
	movq	%rax, 400(%rsp)
	movq	88(%rsp), %rax
	movq	%rax, 408(%rsp)
	testl	%r14d, %r14d
	jne	.L285
	jmp	.L286
.L325:
#APP
# 1140 "../includes/int_double12.0.h" 1
	movapd 0(%rbp),%xmm0
	shufpd $1,%xmm0,%xmm0
	movapd %xmm0,688(%rsp)
	movapd 16(%rbp),%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm1,704(%rsp)
	
# 0 "" 2
#NO_APP
	mov	%r12d, %eax
	movq	%rax, 56(%rsp)
	fildll	56(%rsp)
	fstpl	728(%rsp)
#APP
# 405 "../includes/int_double12.0.h" 1
	movddup 728(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd 0(%r13),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	fldl	448(%rsp)
	fstpl	24(%rsp)
	fldl	456(%rsp)
	fchs
	fstpl	16(%rsp)
	movsd	16(%rsp), %xmm0
	call	log_ru
	movsd	%xmm0, 16(%rsp)
	movsd	24(%rsp), %xmm0
	call	log_rd
	movsd	%xmm0, 8(%rsp)
	fldl	8(%rsp)
	fldl	16(%rsp)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	jbe	.L375
	movapd	%xmm0, %xmm1
	movsd	16(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L375:
	fldl	16(%rsp)
	fchs
	fstpl	472(%rsp)
	movsd	%xmm0, 464(%rsp)
#APP
# 524 "../includes/int_double12.0.h" 1
	movapd 704(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 464(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,448(%rsp)
	
# 0 "" 2
#NO_APP
	fldl	456(%rsp)
	fldl	448(%rsp)
	fxch	%st(1)
	fstl	1416(%rsp)
	fxch	%st(1)
	fstl	1408(%rsp)
	faddp	%st, %st(1)
	fldl	.LC39(%rip)
	fucomip	%st(1), %st
	fstp	%st(0)
	jb	.L376
	movl	$.LC40, %edi
	call	puts
	movl	$1, %edi
	call	exit
.L376:
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd d_pi(%rip),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,536(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 536(%rsp)
	je	.L291
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	movl	$d_pi, %edi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L291:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd d_pi(%rip),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 1408(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,448(%rsp)
# 0 "" 2
#NO_APP
	fldl	456(%rsp)
	fstpl	16(%rsp)
	fldl	448(%rsp)
	fstpl	32(%rsp)
	movsd	32(%rsp), %xmm0
	call	cospi_rd
	movsd	%xmm0, 24(%rsp)
	fldl	16(%rsp)
	fchs
	fstpl	16(%rsp)
	movsd	16(%rsp), %xmm0
	call	cospi_rd
	movsd	%xmm0, 40(%rsp)
	movsd	32(%rsp), %xmm0
	call	sinpi_rd
	movsd	%xmm0, 32(%rsp)
	movsd	16(%rsp), %xmm0
	call	sinpi_rd
	movsd	%xmm0, 16(%rsp)
	fldz
	fldl	24(%rsp)
	fucomip	%st(1), %st
	sbbl	%esi, %esi
	notl	%esi
	andl	$8, %esi
	movl	%esi, %eax
	orl	$4, %eax
	fldl	40(%rsp)
	fucomip	%st(1), %st
	cmovae	%eax, %esi
	movl	%esi, %eax
	orl	$2, %eax
	fldl	32(%rsp)
	fucomip	%st(1), %st
	cmovae	%eax, %esi
	leal	1(%rsi), %eax
	fldl	16(%rsp)
	fucomip	%st(1), %st
	fstp	%st(0)
	cmovae	%eax, %esi
	cmpb	$15, %sil
	ja	.L301
	movzbl	%sil, %eax
	jmp	*(%r15,%rax,8)
	.section	.rodata
	.align 8
	.align 4
.L310:
	.quad	.L302
	.quad	.L301
	.quad	.L303
	.quad	.L304
	.quad	.L305
	.quad	.L301
	.quad	.L301
	.quad	.L301
	.quad	.L301
	.quad	.L301
	.quad	.L301
	.quad	.L306
	.quad	.L307
	.quad	.L308
	.quad	.L301
	.quad	.L309
	.text
.L302:
	fldl	16(%rsp)
	fstpl	336(%rsp)
	fldl	32(%rsp)
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 344(%rsp)
	fldl	24(%rsp)
	fstpl	320(%rsp)
	fldl	40(%rsp)
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 328(%rsp)
	jmp	.L311
.L303:
	fldl	16(%rsp)
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double_neg(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	448(%rsp), %rax
	movq	%rax, 336(%rsp)
	fldl	32(%rsp)
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 344(%rsp)
	fld1
	fchs
	fstpl	320(%rsp)
	fldl	24(%rsp)
	fldl	40(%rsp)
	fxch	%st(1)
	fucomi	%st(1), %st
	fstp	%st(1)
	jb	.L377
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 328(%rsp)
	jmp	.L311
.L377:
	fstp	%st(0)
	fldl	40(%rsp)
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 328(%rsp)
	jmp	.L311
.L304:
	fldl	16(%rsp)
	fstpl	336(%rsp)
	fldl	32(%rsp)
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 344(%rsp)
	fldl	40(%rsp)
	fstpl	320(%rsp)
	fldl	24(%rsp)
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 328(%rsp)
	jmp	.L311
.L305:
	fld1
	fchs
	fstpl	336(%rsp)
	fldl	16(%rsp)
	fldl	32(%rsp)
	fxch	%st(1)
	fucomi	%st(1), %st
	fstp	%st(1)
	jb	.L378
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 344(%rsp)
	jmp	.L316
.L378:
	fstp	%st(0)
	fldl	32(%rsp)
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 344(%rsp)
.L316:
	fldl	24(%rsp)
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double_neg(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	448(%rsp), %rax
	movq	%rax, 320(%rsp)
	fldl	40(%rsp)
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 328(%rsp)
	jmp	.L311
.L306:
	fldl	16(%rsp)
	fldl	32(%rsp)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	movq	16(%rsp), %rax
	cmovae	32(%rsp), %rax
	movq	%rax, 336(%rsp)
	fld1
	fchs
	fstpl	344(%rsp)
	fldl	40(%rsp)
	fstpl	320(%rsp)
	fldl	24(%rsp)
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 328(%rsp)
	jmp	.L311
.L307:
	fldl	32(%rsp)
	fstpl	336(%rsp)
	fldl	16(%rsp)
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 344(%rsp)
	fldl	24(%rsp)
	fstpl	320(%rsp)
	fldl	40(%rsp)
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 328(%rsp)
	jmp	.L311
.L308:
	fld1
	fchs
	fstpl	328(%rsp)
	fldl	40(%rsp)
	fldl	24(%rsp)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	movq	40(%rsp), %rax
	cmovae	24(%rsp), %rax
	movq	%rax, 320(%rsp)
	fldl	32(%rsp)
	fstpl	336(%rsp)
	fldl	16(%rsp)
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 344(%rsp)
	jmp	.L311
.L309:
	fldl	32(%rsp)
	fstpl	336(%rsp)
	fldl	16(%rsp)
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 344(%rsp)
	fldl	40(%rsp)
	fstpl	320(%rsp)
	fldl	24(%rsp)
	fstpl	536(%rsp)
#APP
# 418 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd delta_int_double(%rip),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 328(%rsp)
	jmp	.L311
.L301:
	movzbl	%sil, %esi
	movl	$.LC41, %edi
	movl	$0, %eax
	call	printf
	movsd	40(%rsp), %xmm1
	movsd	24(%rsp), %xmm0
	movl	$.LC42, %edi
	movl	$2, %eax
	call	printf
	movsd	16(%rsp), %xmm1
	movsd	32(%rsp), %xmm0
	movl	$.LC43, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L311:
#APP
# 524 "../includes/int_double12.0.h" 1
	movapd 688(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 464(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,448(%rsp)
	
# 0 "" 2
#NO_APP
	fldl	448(%rsp)
	fstpl	24(%rsp)
	fldl	456(%rsp)
	fchs
	fstpl	16(%rsp)
	movsd	16(%rsp), %xmm0
	call	exp_ru
	movsd	%xmm0, 16(%rsp)
	movsd	24(%rsp), %xmm0
	call	exp_rd
	movsd	%xmm0, 8(%rsp)
	fldl	8(%rsp)
	fldl	16(%rsp)
	fxch	%st(1)
	fucomip	%st(1), %st
	fstp	%st(0)
	jbe	.L379
	movapd	%xmm0, %xmm1
	movsd	16(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L379:
	fldl	16(%rsp)
	fchs
	fstpl	1432(%rsp)
	movsd	%xmm0, 1424(%rsp)
#APP
# 524 "../includes/int_double12.0.h" 1
	movapd 1424(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 336(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rdx
	movq	448(%rsp), %rax
#APP
# 524 "../includes/int_double12.0.h" 1
	movapd 1424(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 320(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rcx
	movq	%rcx, 664(%rsp)
	movq	448(%rsp), %rcx
	movq	%rcx, 656(%rsp)
	movq	%rdx, 680(%rsp)
	movq	%rax, 672(%rsp)
#APP
# 1048 "../includes/int_double12.0.h" 1
	movapd (%rbx),%xmm0
	addpd 656(%rsp),%xmm0
	movapd %xmm0,1184(%rsp)
	movapd 16(%rbx),%xmm1
	addpd 672(%rsp),%xmm1
	movapd %xmm1,1200(%rsp)
	
# 0 "" 2
#NO_APP
	movq	1184(%rsp), %rax
	movq	%rax, (%rbx)
	movq	1192(%rsp), %rax
	movq	%rax, 8(%rbx)
	movq	1200(%rsp), %rax
	movq	%rax, 16(%rbx)
	movq	1208(%rsp), %rax
	movq	%rax, 24(%rbx)
	addl	$1, %r12d
	cmpl	%r14d, %r12d
	jb	.L325
.L286:
	fld1
	fstpl	808(%rsp)
#APP
# 1187 "../includes/int_double12.0.h" 1
	movddup 808(%rsp),%xmm0
	movapd 0(%rbp),%xmm1
	addsubpd %xmm0,%xmm1
	movapd %xmm1,768(%rsp)
	movapd 16(%rbp),%xmm1
	movapd %xmm1,784(%rsp)
	
# 0 "" 2
# 961 "../includes/int_double12.0.h" 1
	movapd 784(%rsp),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,448(%rsp)
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 1448(%rsp)
	movq	448(%rsp), %rax
	movq	%rax, 1440(%rsp)
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd 768(%rsp),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,448(%rsp)
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 1464(%rsp)
	movq	448(%rsp), %rax
	movq	%rax, 1456(%rsp)
#APP
# 392 "../includes/int_double12.0.h" 1
	movapd 1456(%rsp),%XMM0
	addpd 1440(%rsp),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 1544(%rsp)
	movq	448(%rsp), %rax
	movq	%rax, 1536(%rsp)
#APP
# 1033 "../includes/int_double12.0.h" 1
	movapd 768(%rsp),%xmm0
	movapd %xmm0,1472(%rsp)
	movapd 784(%rsp),%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm1,1488(%rsp)
	
# 0 "" 2
# 1318 "../includes/int_double12.0.h" 1
	movapd 400(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1488(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	shufpd $1,%xmm7,%xmm7
	movapd 384(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1472(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm7,%xmm0
	movapd %xmm0,1504(%rsp)
	movapd 400(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1472(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	movapd 384(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1488(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm0,%xmm7
	movapd %xmm7,1520(%rsp)
# 0 "" 2
# 1003 "../includes/int_double12.0.h" 1
	movapd 1536(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,536(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 536(%rsp)
	je	.L326
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	1536(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L326:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 1536(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 1520(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,448(%rsp)
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rsi
	movq	448(%rsp), %rdx
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd 1536(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,536(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 536(%rsp)
	je	.L327
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	1536(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L327:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 1536(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 1504(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,448(%rsp)
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 744(%rsp)
	movq	448(%rsp), %rax
	movq	%rax, 736(%rsp)
	movq	%rsi, 760(%rsp)
	movq	%rdx, 752(%rsp)
#APP
# 1048 "../includes/int_double12.0.h" 1
	movapd (%rbx),%xmm0
	addpd 736(%rsp),%xmm0
	movapd %xmm0,1216(%rsp)
	movapd 16(%rbx),%xmm1
	addpd 752(%rsp),%xmm1
	movapd %xmm1,1232(%rsp)
	
# 0 "" 2
#NO_APP
	movq	1216(%rsp), %rax
	movq	%rax, (%rbx)
	movq	1224(%rsp), %rax
	movq	%rax, 8(%rbx)
	movq	1232(%rsp), %rax
	movq	%rax, 16(%rbx)
	movq	1240(%rsp), %rax
	movq	%rax, 24(%rbx)
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd 496(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,536(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 536(%rsp)
	je	.L328
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	496(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L328:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 496(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 400(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,448(%rsp)
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rsi
	movq	448(%rsp), %rdx
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd 496(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,536(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 536(%rsp)
	je	.L329
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	496(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L329:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 496(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 384(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,448(%rsp)
# 0 "" 2
#NO_APP
	movq	%rsi, 408(%rsp)
	movq	%rdx, 400(%rsp)
	movq	456(%rsp), %rax
	movq	%rax, 392(%rsp)
	movq	448(%rsp), %rax
	movq	%rax, 384(%rsp)
	flds	.LC19(%rip)
	fstpl	856(%rsp)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup 856(%rsp),%xmm0
	movapd 400(%rsp),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,448(%rsp)
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rdx
	movq	448(%rsp), %rax
	fldl	856(%rsp)
	fldz
	fxch	%st(1)
	fucomi	%st(1), %st
	fstp	%st(1)
	jbe	.L380
	fstp	%st(0)
#APP
# 640 "../includes/int_double12.0.h" 1
	movddup 856(%rsp),%xmm0
	movapd 384(%rsp),%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,448(%rsp)
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rsi
	movq	448(%rsp), %rcx
	jmp	.L332
.L380:
	fldz
	fucomip	%st(1), %st
	jbe	.L381
	fchs
	fstpl	536(%rsp)
#APP
# 653 "../includes/int_double12.0.h" 1
	movddup 536(%rsp),%xmm0
	movapd 384(%rsp),%xmm1
	shufpd $1,%xmm1,%xmm1
	divpd %xmm0,%xmm1
	movapd %xmm1,448(%rsp)
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rsi
	movq	448(%rsp), %rcx
	jmp	.L332
.L381:
	fstp	%st(0)
	movl	$.LC49, %edi
	call	puts
	movl	$1, %edi
	call	exit
.L332:
	movq	%rsi, 824(%rsp)
	movq	%rcx, 816(%rsp)
	movq	%rdx, 840(%rsp)
	movq	%rax, 832(%rsp)
#APP
# 1048 "../includes/int_double12.0.h" 1
	movapd (%rbx),%xmm0
	addpd 816(%rsp),%xmm0
	movapd %xmm0,1248(%rsp)
	movapd 16(%rbx),%xmm1
	addpd 832(%rsp),%xmm1
	movapd %xmm1,1264(%rsp)
	
# 0 "" 2
#NO_APP
	movq	1248(%rsp), %rax
	movq	%rax, (%rbx)
	movq	1256(%rsp), %rax
	movq	%rax, 8(%rbx)
	movq	1264(%rsp), %rax
	movq	%rax, 16(%rbx)
	movq	1272(%rsp), %rax
	movq	%rax, 24(%rbx)
#APP
# 1318 "../includes/int_double12.0.h" 1
	movapd 16(%rbp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 400(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	shufpd $1,%xmm7,%xmm7
	movapd 0(%rbp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 384(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm7,%xmm0
	movapd %xmm0,864(%rsp)
	movapd 16(%rbp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 384(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	movapd 0(%rbp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 400(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm0,%xmm7
	movapd %xmm7,880(%rsp)
# 0 "" 2
# 1003 "../includes/int_double12.0.h" 1
	movapd 496(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,536(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 536(%rsp)
	je	.L335
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	496(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L335:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 496(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 880(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,448(%rsp)
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rsi
	movq	448(%rsp), %rdx
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd 496(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,536(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 536(%rsp)
	je	.L336
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	496(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L336:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 496(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 864(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,448(%rsp)
# 0 "" 2
#NO_APP
	movq	%rsi, 120(%rsp)
	movq	%rdx, 112(%rsp)
	movq	456(%rsp), %rax
	movq	%rax, 104(%rsp)
	movq	448(%rsp), %rax
	movq	%rax, 96(%rsp)
	fld1
	fstpl	904(%rsp)
#APP
# 1091 "../includes/int_double12.0.h" 1
	movddup 904(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	addpd 416(%rsp),%xmm0
	movapd %xmm0,1280(%rsp)
	movapd 432(%rsp),%xmm1
	movapd %xmm1,1296(%rsp)
	
# 0 "" 2
#NO_APP
	movq	1280(%rsp), %rax
	movq	%rax, 416(%rsp)
	movq	1288(%rsp), %rax
	movq	%rax, 424(%rsp)
	movq	1296(%rsp), %rax
	movq	%rax, 432(%rsp)
	movq	1304(%rsp), %rax
	movq	%rax, 440(%rsp)
#APP
# 1318 "../includes/int_double12.0.h" 1
	movapd 112(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 432(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	shufpd $1,%xmm7,%xmm7
	movapd 96(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 416(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm7,%xmm0
	movapd %xmm0,912(%rsp)
	movapd 112(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 416(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	movapd 96(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 432(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm0,%xmm7
	movapd %xmm7,928(%rsp)
# 0 "" 2
# 1003 "../includes/int_double12.0.h" 1
	movapd 496(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,536(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 536(%rsp)
	je	.L382
	jmp	.L337
.L343:
	fxch	%st(1)
	fstl	904(%rsp)
#APP
# 1091 "../includes/int_double12.0.h" 1
	movddup 904(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	addpd 416(%rsp),%xmm0
	movapd %xmm0,1280(%rsp)
	movapd 432(%rsp),%xmm1
	movapd %xmm1,1296(%rsp)
	
# 0 "" 2
#NO_APP
	movq	1280(%rsp), %rax
	movq	%rax, 416(%rsp)
	movq	1288(%rsp), %rax
	movq	%rax, 424(%rsp)
	movq	1296(%rsp), %rax
	movq	%rax, 432(%rsp)
	movq	1304(%rsp), %rax
	movq	%rax, 440(%rsp)
#APP
# 1318 "../includes/int_double12.0.h" 1
	movapd 16(%rsi),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 432(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	shufpd $1,%xmm7,%xmm7
	movapd (%rsi),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 416(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm7,%xmm0
	movapd %xmm0,912(%rsp)
	movapd 16(%rsi),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 416(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	movapd (%rsi),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 432(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm0,%xmm7
	movapd %xmm7,928(%rsp)
# 0 "" 2
# 1003 "../includes/int_double12.0.h" 1
	movapd 496(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,536(%rsp)
# 0 "" 2
#NO_APP
	addl	$1, %edi
	addq	$32, %rdx
	cmpl	$0, 536(%rsp)
	je	.L339
	fstp	%st(0)
	fstp	%st(0)
.L337:
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	496(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L382:
	leaq	152(%rsp), %rdx
	movl	$1, %edi
	fld1
	leaq	96(%rsp), %r9
	fld	%st(0)
.L339:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 496(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 928(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,448(%rsp)
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rbp
	movq	448(%rsp), %rsi
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd 496(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,536(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 536(%rsp)
	je	.L340
	fstp	%st(0)
	fstp	%st(0)
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	496(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L340:
	fxch	%st(1)
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 496(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 912(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,448(%rsp)
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rcx
	movq	448(%rsp), %rax
	movq	%rbp, (%rdx)
	movq	%rsi, -8(%rdx)
	movq	%rcx, -16(%rdx)
	movq	%rax, -24(%rdx)
	fstl	952(%rsp)
#APP
# 1091 "../includes/int_double12.0.h" 1
	movddup 952(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	addpd 416(%rsp),%xmm0
	movapd %xmm0,1312(%rsp)
	movapd 432(%rsp),%xmm1
	movapd %xmm1,1328(%rsp)
	
# 0 "" 2
#NO_APP
	movq	1312(%rsp), %rax
	movq	%rax, 416(%rsp)
	movq	1320(%rsp), %rax
	movq	%rax, 424(%rsp)
	movq	1328(%rsp), %rax
	movq	%rax, 432(%rsp)
	movq	1336(%rsp), %rax
	movq	%rax, 440(%rsp)
	mov	%edi, %esi
	salq	$5, %rsi
	leaq	(%r9,%rsi), %rsi
#APP
# 1318 "../includes/int_double12.0.h" 1
	movapd 16(%rsi),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 432(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	shufpd $1,%xmm7,%xmm7
	movapd (%rsi),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 416(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm7,%xmm0
	movapd %xmm0,960(%rsp)
	movapd 16(%rsi),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 416(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,%xmm7
	movapd (%rsi),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 432(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	addpd %xmm0,%xmm7
	movapd %xmm7,976(%rsp)
# 0 "" 2
# 1003 "../includes/int_double12.0.h" 1
	movapd 496(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,536(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 536(%rsp)
	je	.L341
	fstp	%st(0)
	fstp	%st(0)
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	496(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L341:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 496(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 976(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,448(%rsp)
# 0 "" 2
#NO_APP
	movq	456(%rsp), %r8
	movq	448(%rsp), %rbp
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd 496(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,536(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 536(%rsp)
	je	.L342
	fstp	%st(0)
	fstp	%st(0)
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	496(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L342:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 496(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 960(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,448(%rsp)
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rcx
	movq	448(%rsp), %rax
	movq	%r8, (%rdx)
	movq	%rbp, -8(%rdx)
	movq	%rcx, -16(%rdx)
	movq	%rax, -24(%rdx)
	cmpl	$6, %edi
	jne	.L343
	fstp	%st(0)
	fstp	%st(0)
	movl	$0, %eax
	leaq	96(%rsp), %rsi
.L344:
	mov	%eax, %edx
	movq	%rdx, %rcx
	salq	$4, %rcx
	addq	$h_bernoulli, %rcx
	leaq	1(%rdx,%rdx), %rdi
	salq	$4, %rdi
#APP
# 524 "../includes/int_double12.0.h" 1
	movapd (%rcx),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd (%rsi,%rdi),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rdi
	movq	448(%rsp), %rbp
	salq	$5, %rdx
#APP
# 524 "../includes/int_double12.0.h" 1
	movapd (%rcx),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd (%rsi,%rdx),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	%rbp, 368(%rsp)
	movq	%rdi, 376(%rsp)
	movq	448(%rsp), %rdx
	movq	%rdx, 352(%rsp)
	movq	456(%rsp), %rdx
	movq	%rdx, 360(%rsp)
#APP
# 1048 "../includes/int_double12.0.h" 1
	movapd (%rbx),%xmm0
	addpd 352(%rsp),%xmm0
	movapd %xmm0,1344(%rsp)
	movapd 16(%rbx),%xmm1
	addpd 368(%rsp),%xmm1
	movapd %xmm1,1360(%rsp)
	
# 0 "" 2
#NO_APP
	movq	1344(%rsp), %rdx
	movq	%rdx, (%rbx)
	movq	1352(%rsp), %rdx
	movq	%rdx, 8(%rbx)
	movq	1360(%rsp), %rdx
	movq	%rdx, 16(%rbx)
	movq	1368(%rsp), %rdx
	movq	%rdx, 24(%rbx)
	addl	$1, %eax
	cmpl	$7, %eax
	jne	.L344
	fld1
	fstl	1032(%rsp)
	flds	.LC87(%rip)
	fstl	1064(%rsp)
	fxch	%st(1)
#APP
# 405 "../includes/int_double12.0.h" 1
	movddup 1064(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd 480(%rsp),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 1048(%rsp)
	movq	448(%rsp), %rax
	movq	%rax, 1040(%rsp)
#APP
# 467 "../includes/int_double12.0.h" 1
	movddup 1032(%rsp),%XMM0
	movapd 1040(%rsp),%XMM1
	addsubpd %XMM0,%XMM1
	movapd %XMM1,448(%rsp)
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 1016(%rsp)
	movq	448(%rsp), %rax
	movq	%rax, 1008(%rsp)
	fstpl	1096(%rsp)
	fstpl	1128(%rsp)
	movq	48(%rsp), %rax
#APP
# 405 "../includes/int_double12.0.h" 1
	movddup 1128(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd (%rax),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 1112(%rsp)
	movq	448(%rsp), %rax
	movq	%rax, 1104(%rsp)
#APP
# 467 "../includes/int_double12.0.h" 1
	movddup 1096(%rsp),%XMM0
	movapd 1104(%rsp),%XMM1
	addsubpd %XMM0,%XMM1
	movapd %XMM1,448(%rsp)
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 1080(%rsp)
	movq	448(%rsp), %rax
	movq	%rax, 1072(%rsp)
#APP
# 524 "../includes/int_double12.0.h" 1
	movapd 1008(%rsp),%xmm0
	xorpd	d_zero(%rip),%xmm0
	movapd 1072(%rsp),%xmm1
	movapd %xmm1,%xmm2
	shufpd	$1,%xmm2,%xmm2
	movapd %xmm0,%xmm3
	cmpltpd d_neg_zero(%rip),%xmm3
	xorpd	d_neg_zero(%rip),%xmm2
	movapd %xmm3,%xmm4
	shufpd $1,%xmm4,%xmm4
	movapd %xmm3,%xmm5
	andpd %xmm2,%xmm5
	andnpd %xmm1,%xmm3
	movapd %xmm4,%xmm6
	andnpd %xmm1,%xmm6
	andpd %xmm4,%xmm2
	orpd %xmm5,%xmm3
	movapd %xmm0,%xmm4
	shufpd $1,%xmm4,%xmm4
	orpd %xmm6,%xmm2
	mulpd %xmm3,%xmm0
	mulpd %xmm4,%xmm2
	minpd %xmm2,%xmm0
	movapd %xmm0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 1000(%rsp)
	movq	448(%rsp), %rax
	movq	%rax, 992(%rsp)
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd 368(%rsp),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,448(%rsp)
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 1400(%rsp)
	movq	448(%rsp), %rax
	movq	%rax, 1392(%rsp)
#APP
# 961 "../includes/int_double12.0.h" 1
	movapd 352(%rsp),%xmm0
	movapd %xmm0,%xmm1
	movapd %xmm0,%xmm2
	shufpd $1,%xmm1,%xmm1
	minsd %xmm1,%xmm0
	maxsd %xmm2,%xmm1
	maxsd d_zero_zero(%rip),%xmm1
	unpcklpd %xmm0,%xmm1
	movapd %xmm1,%xmm0
	xorpd d_zero(%rip),%xmm1
	mulpd %xmm1,%xmm0
	movapd %xmm0,448(%rsp)
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 1384(%rsp)
	movq	448(%rsp), %rax
	movq	%rax, 1376(%rsp)
#APP
# 392 "../includes/int_double12.0.h" 1
	movapd 1376(%rsp),%XMM0
	addpd 1392(%rsp),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	448(%rsp), %rax
	movq	%rax, 560(%rsp)
	movq	456(%rsp), %rax
	movq	%rax, 552(%rsp)
#APP
# 979 "../includes/int_double12.0.h" 1
	fldl 560(%rsp)
	fsqrt
	fstpl 544(%rsp)
	fldl 552(%rsp)
	fchs
	fsqrt
	fstpl 536(%rsp)
	
# 0 "" 2
#NO_APP
	fldl	536(%rsp)
	fldl	544(%rsp)
	fucomi	%st(1), %st
	jbe	.L383
	fstpl	8(%rsp)
	movsd	8(%rsp), %xmm1
	fstpl	8(%rsp)
	movsd	8(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L383:
	fxch	%st(1)
	fchs
	fstpl	456(%rsp)
	fstpl	448(%rsp)
#APP
# 392 "../includes/int_double12.0.h" 1
	movapd 448(%rsp),%XMM0
	addpd delta_blow(%rip),%XMM0
	movapd %XMM0,464(%rsp)
	
# 0 "" 2
#NO_APP
	movq	472(%rsp), %rax
	movq	%rax, 1144(%rsp)
	movq	464(%rsp), %rax
	movq	%rax, 1136(%rsp)
#APP
# 1003 "../includes/int_double12.0.h" 1
	movapd 992(%rsp),%xmm0
	cmplepd d_zero_zero(%rip),%xmm0
	movmskpd %xmm0,%eax
	movl %eax,%ecx
	shr %eax
	and %ecx,%eax
	movl %eax,536(%rsp)
# 0 "" 2
#NO_APP
	cmpl	$0, 536(%rsp)
	je	.L347
	movl	$.LC37, %edi
	call	puts
	movl	$32, %edi
	call	putchar
	leaq	992(%rsp), %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$1, %edi
	call	exit
.L347:
#APP
# 619 "../includes/int_double12.0.h" 1
	movapd 992(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	movapd %xmm0,%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm0,%xmm2
	cmpnlepd d_zero_zero(%rip),%xmm2
	movapd 1136(%rsp),%xmm3
	movapd %xmm3,%xmm4
	shufpd $1,%xmm3,%xmm3
	xorpd d_neg_zero(%rip),%xmm3
	andpd %xmm2,%xmm4
	andnpd %xmm3,%xmm2
	orpd %xmm4,%xmm2
	movapd %xmm2,%xmm3
	divpd %xmm1,%xmm2
	divpd %xmm0,%xmm3
	minpd %xmm2,%xmm3
	movapd %xmm3,448(%rsp)
# 0 "" 2
#NO_APP
	fldl	448(%rsp)
	movq	456(%rsp), %rax
	movq	%rax, 520(%rsp)
	fstl	512(%rsp)
	movq	%rax, 8(%rsp)
	fldl	8(%rsp)
	fucomi	%st(1), %st
	jb	.L384
	fstp	%st(0)
	fstpl	520(%rsp)
	jmp	.L350
.L384:
	fstp	%st(1)
	fstpl	512(%rsp)
.L350:
#APP
# 392 "../includes/int_double12.0.h" 1
	movapd (%rbx),%XMM0
	addpd 512(%rsp),%XMM0
	movapd %XMM0,464(%rsp)
	
# 0 "" 2
#NO_APP
	movq	472(%rsp), %rax
	movq	%rax, 8(%rbx)
	movq	464(%rsp), %rax
	movq	%rax, (%rbx)
#APP
# 392 "../includes/int_double12.0.h" 1
	movapd 16(%rbx),%XMM0
	addpd 512(%rsp),%XMM0
	movapd %XMM0,464(%rsp)
	
# 0 "" 2
#NO_APP
	movq	472(%rsp), %rax
	movq	%rax, 24(%rbx)
	movq	464(%rsp), %rax
	movq	%rax, 16(%rbx)
	movq	%rbx, %rax
	addq	$1560, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
.L374:
	.cfi_restore_state
	fstp	%st(0)
	flds	.LC86(%rip)
	fstpl	568(%rsp)
#APP
# 405 "../includes/int_double12.0.h" 1
	movddup 568(%rsp),%XMM0
	xorpd d_zero(%rip),%XMM0
	addpd 0(%r13),%XMM0
	movapd %XMM0,448(%rsp)
	
# 0 "" 2
#NO_APP
	movq	456(%rsp), %rax
	movq	%rax, 504(%rsp)
	movq	448(%rsp), %rax
	movq	%rax, 496(%rsp)
	fld1
	fstpl	616(%rsp)
#APP
# 1140 "../includes/int_double12.0.h" 1
	movapd 0(%rbp),%xmm0
	shufpd $1,%xmm0,%xmm0
	movapd %xmm0,624(%rsp)
	movapd 16(%rbp),%xmm1
	shufpd $1,%xmm1,%xmm1
	movapd %xmm1,640(%rsp)
	
# 0 "" 2
# 1091 "../includes/int_double12.0.h" 1
	movddup 616(%rsp),%xmm0
	xorpd d_zero(%rip),%xmm0
	addpd 624(%rsp),%xmm0
	movapd %xmm0,576(%rsp)
	movapd 640(%rsp),%xmm1
	movapd %xmm1,592(%rsp)
	
# 0 "" 2
#NO_APP
	leaq	64(%rsp), %rdi
	leaq	576(%rsp), %rdx
	leaq	496(%rsp), %rsi
	call	_Z3powRK10int_doubleRK11int_complex
	movq	64(%rsp), %rax
	movq	%rax, 384(%rsp)
	movq	72(%rsp), %rax
	movq	%rax, 392(%rsp)
	movq	80(%rsp), %rax
	movq	%rax, 400(%rsp)
	movq	88(%rsp), %rax
	movq	%rax, 408(%rsp)
	movl	$50, %r14d
.L285:
	movl	$0, %r12d
	movl	$.L310, %r15d
	jmp	.L325
	.cfi_endproc
.LFE1163:
	.size	_Z7hurwitzRK11int_complexRK10int_double, .-_Z7hurwitzRK11int_complexRK10int_double
.globl _i_pi
	.data
	.align 4
	.type	_i_pi, @object
	.size	_i_pi, 8
_i_pi:
	.long	1413754136
	.long	1074340347
.globl _i_pi2
	.align 4
	.type	_i_pi2, @object
	.size	_i_pi2, 8
_i_pi2:
	.long	1413754137
	.long	1074340347
.globl _d_pi
	.align 8
	.type	_d_pi, @object
	.size	_d_pi, 8
_d_pi:
	.quad	_i_pi
.globl _d_pi2
	.align 8
	.type	_d_pi2, @object
	.size	_d_pi2, 8
_d_pi2:
	.quad	_i_pi2
.globl d_pi
	.bss
	.align 16
	.type	d_pi, @object
	.size	d_pi, 16
d_pi:
	.zero	16
.globl _i_gamma
	.data
	.align 4
	.type	_i_gamma, @object
	.size	_i_gamma, 8
_i_gamma:
	.long	-59787752
	.long	1071806604
.globl _i_gamma2
	.align 4
	.type	_i_gamma2, @object
	.size	_i_gamma2, 8
_i_gamma2:
	.long	-59787751
	.long	1071806604
.globl _d_gamma
	.align 8
	.type	_d_gamma, @object
	.size	_d_gamma, 8
_d_gamma:
	.quad	_i_gamma
.globl _d_gamma2
	.align 8
	.type	_d_gamma2, @object
	.size	_d_gamma2, 8
_d_gamma2:
	.quad	_i_gamma2
.globl d_gamma
	.bss
	.align 16
	.type	d_gamma, @object
	.size	d_gamma, 16
d_gamma:
	.zero	16
.globl _i_pi_2
	.data
	.align 4
	.type	_i_pi_2, @object
	.size	_i_pi_2, 8
_i_pi_2:
	.long	1413754136
	.long	1073291771
.globl _i_pi2_2
	.align 4
	.type	_i_pi2_2, @object
	.size	_i_pi2_2, 8
_i_pi2_2:
	.long	1413754137
	.long	1073291771
.globl _d_pi_2
	.align 8
	.type	_d_pi_2, @object
	.size	_d_pi_2, 8
_d_pi_2:
	.quad	_i_pi_2
.globl _d_pi2_2
	.align 8
	.type	_d_pi2_2, @object
	.size	_d_pi2_2, 8
_d_pi2_2:
	.quad	_i_pi2_2
.globl d_pi_2
	.bss
	.align 16
	.type	d_pi_2, @object
	.size	d_pi_2, 16
d_pi_2:
	.zero	16
.globl _i_2pi_2
	.data
	.align 4
	.type	_i_2pi_2, @object
	.size	_i_2pi_2, 8
_i_2pi_2:
	.long	1413754136
	.long	1075388923
.globl _i_2pi2_2
	.align 4
	.type	_i_2pi2_2, @object
	.size	_i_2pi2_2, 8
_i_2pi2_2:
	.long	1413754137
	.long	1075388923
.globl _d_2pi_2
	.align 8
	.type	_d_2pi_2, @object
	.size	_d_2pi_2, 8
_d_2pi_2:
	.quad	_i_2pi_2
.globl _d_2pi2_2
	.align 8
	.type	_d_2pi2_2, @object
	.size	_d_2pi2_2, 8
_d_2pi2_2:
	.quad	_i_2pi2_2
.globl d_two_pi
	.bss
	.align 16
	.type	d_two_pi, @object
	.size	d_two_pi, 16
d_two_pi:
	.zero	16
.globl d_ln_pi
	.align 16
	.type	d_ln_pi, @object
	.size	d_ln_pi, 16
d_ln_pi:
	.zero	16
.globl d_ln_two_pi
	.align 16
	.type	d_ln_two_pi, @object
	.size	d_ln_two_pi, 16
d_ln_two_pi:
	.zero	16
.globl __nze
	.data
	.align 4
	.type	__nze, @object
	.size	__nze, 8
__nze:
	.long	0
	.long	-2147483648
.globl _nze
	.bss
	.align 8
	.type	_nze, @object
	.size	_nze, 8
_nze:
	.zero	8
.globl d_zero_zero
	.align 16
	.type	d_zero_zero, @object
	.size	d_zero_zero, 16
d_zero_zero:
	.zero	16
.globl d_zero
	.align 16
	.type	d_zero, @object
	.size	d_zero, 16
d_zero:
	.zero	16
.globl d_neg_zero
	.align 16
	.type	d_neg_zero, @object
	.size	d_neg_zero, 16
d_neg_zero:
	.zero	16
.globl d_neg_neg_zero
	.align 16
	.type	d_neg_neg_zero, @object
	.size	d_neg_neg_zero, 16
d_neg_neg_zero:
	.zero	16
.globl d_one
	.align 16
	.type	d_one, @object
	.size	d_one, 16
d_one:
	.zero	16
.globl c_zero
	.align 32
	.type	c_zero, @object
	.size	c_zero, 32
c_zero:
	.zero	32
.globl c_one
	.align 32
	.type	c_one, @object
	.size	c_one, 32
c_one:
	.zero	32
.globl c_half
	.align 32
	.type	c_half, @object
	.size	c_half, 32
c_half:
	.zero	32
.globl delta_int_double
	.align 16
	.type	delta_int_double, @object
	.size	delta_int_double, 16
delta_int_double:
	.zero	16
.globl delta_int_double_neg
	.align 16
	.type	delta_int_double_neg, @object
	.size	delta_int_double_neg, 16
delta_int_double_neg:
	.zero	16
.globl delta_blow
	.align 16
	.type	delta_blow, @object
	.size	delta_blow, 16
delta_blow:
	.zero	16
.globl old_cw
	.align 2
	.type	old_cw, @object
	.size	old_cw, 2
old_cw:
	.zero	2
.globl new_cw
	.align 2
	.type	new_cw, @object
	.size	new_cw, 2
new_cw:
	.zero	2
.globl old__SSE_cw
	.align 4
	.type	old__SSE_cw, @object
	.size	old__SSE_cw, 4
old__SSE_cw:
	.zero	4
.globl new__SSE_cw
	.align 4
	.type	new__SSE_cw, @object
	.size	new__SSE_cw, 4
new__SSE_cw:
	.zero	4
.globl bernoulli
	.align 32
	.type	bernoulli, @object
	.size	bernoulli, 112
bernoulli:
	.zero	112
.globl h_bernoulli
	.align 32
	.type	h_bernoulli, @object
	.size	h_bernoulli, 112
h_bernoulli:
	.zero	112
.globl h_bernoulli_initialised
	.type	h_bernoulli_initialised, @object
	.size	h_bernoulli_initialised, 1
h_bernoulli_initialised:
	.zero	1
.globl ln_gamma_err
	.align 32
	.type	ln_gamma_err, @object
	.size	ln_gamma_err, 32
ln_gamma_err:
	.zero	32
.globl ln_gamma_err_2
	.align 32
	.type	ln_gamma_err_2, @object
	.size	ln_gamma_err_2, 32
ln_gamma_err_2:
	.zero	32
.globl re_ln_gamma_err
	.align 16
	.type	re_ln_gamma_err, @object
	.size	re_ln_gamma_err, 16
re_ln_gamma_err:
	.zero	16
	.local	_ZStL8__ioinit
	.comm	_ZStL8__ioinit,1,1
	.local	_ZL6d_half
	.comm	_ZL6d_half,16,16
	.weakref	_ZL20__gthrw_pthread_oncePiPFvvE,pthread_once
	.weakref	_ZL27__gthrw_pthread_getspecificj,pthread_getspecific
	.weakref	_ZL27__gthrw_pthread_setspecificjPKv,pthread_setspecific
	.weakref	_ZL22__gthrw_pthread_createPmPK14pthread_attr_tPFPvS3_ES3_,pthread_create
	.weakref	_ZL20__gthrw_pthread_joinmPPv,pthread_join
	.weakref	_ZL21__gthrw_pthread_equalmm,pthread_equal
	.weakref	_ZL20__gthrw_pthread_selfv,pthread_self
	.weakref	_ZL22__gthrw_pthread_detachm,pthread_detach
	.weakref	_ZL22__gthrw_pthread_cancelm,pthread_cancel
	.weakref	_ZL19__gthrw_sched_yieldv,sched_yield
	.weakref	_ZL26__gthrw_pthread_mutex_lockP15pthread_mutex_t,pthread_mutex_lock
	.weakref	_ZL29__gthrw_pthread_mutex_trylockP15pthread_mutex_t,pthread_mutex_trylock
	.weakref	_ZL31__gthrw_pthread_mutex_timedlockP15pthread_mutex_tPK8timespec,pthread_mutex_timedlock
	.weakref	_ZL28__gthrw_pthread_mutex_unlockP15pthread_mutex_t,pthread_mutex_unlock
	.weakref	_ZL26__gthrw_pthread_mutex_initP15pthread_mutex_tPK19pthread_mutexattr_t,pthread_mutex_init
	.weakref	_ZL29__gthrw_pthread_mutex_destroyP15pthread_mutex_t,pthread_mutex_destroy
	.weakref	_ZL30__gthrw_pthread_cond_broadcastP14pthread_cond_t,pthread_cond_broadcast
	.weakref	_ZL27__gthrw_pthread_cond_signalP14pthread_cond_t,pthread_cond_signal
	.weakref	_ZL25__gthrw_pthread_cond_waitP14pthread_cond_tP15pthread_mutex_t,pthread_cond_wait
	.weakref	_ZL30__gthrw_pthread_cond_timedwaitP14pthread_cond_tP15pthread_mutex_tPK8timespec,pthread_cond_timedwait
	.weakref	_ZL28__gthrw_pthread_cond_destroyP14pthread_cond_t,pthread_cond_destroy
	.weakref	_ZL26__gthrw_pthread_key_createPjPFvPvE,pthread_key_create
	.weakref	_ZL26__gthrw_pthread_key_deletej,pthread_key_delete
	.weakref	_ZL30__gthrw_pthread_mutexattr_initP19pthread_mutexattr_t,pthread_mutexattr_init
	.weakref	_ZL33__gthrw_pthread_mutexattr_settypeP19pthread_mutexattr_ti,pthread_mutexattr_settype
	.weakref	_ZL33__gthrw_pthread_mutexattr_destroyP19pthread_mutexattr_t,pthread_mutexattr_destroy
	.section	.rodata.cst4,"aM",@progbits,4
	.align 4
.LC3:
	.long	1056964608
	.align 4
.LC8:
	.long	3204448256
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC9:
	.long	0
	.long	1048576
	.align 8
.LC10:
	.long	0
	.long	-2146435072
	.align 8
.LC11:
	.long	1851457824
	.long	1024132487
	.align 8
.LC12:
	.long	1851457824
	.long	-1123351161
	.align 8
.LC13:
	.long	3256091556
	.long	1027038407
	.align 8
.LC14:
	.long	3256091556
	.long	-1120445241
	.section	.rodata.cst4
	.align 4
.LC18:
	.long	1086324736
	.align 4
.LC19:
	.long	1073741824
	.align 4
.LC20:
	.long	1106247680
	.align 4
.LC21:
	.long	1103101952
	.align 4
.LC22:
	.long	1109917696
	.align 4
.LC23:
	.long	1144258560
	.align 4
.LC24:
	.long	1193115648
	.align 4
.LC25:
	.long	1115947008
	.align 4
.LC26:
	.long	1084227584
	.align 4
.LC27:
	.long	3231711232
	.align 4
.LC28:
	.long	1247640576
	.align 4
.LC29:
	.long	1160421376
	.align 4
.LC30:
	.long	3291267072
	.align 4
.LC31:
	.long	1143783424
	.align 4
.LC32:
	.long	1306814432
	.align 4
.LC33:
	.long	1088421888
	.align 4
.LC34:
	.long	3235905536
	.align 4
.LC35:
	.long	1127612416
	.align 4
.LC36:
	.long	1097859072
	.align 4
.LC38:
	.long	1602224128
	.section	.rodata.cst8
	.align 8
.LC39:
	.long	1413754136
	.long	-1074191877
	.section	.rodata.cst4
	.align 4
.LC50:
	.long	3196059648
	.align 4
.LC51:
	.long	1048576000
	.align 4
.LC52:
	.long	1092616192
	.align 4
.LC54:
	.long	2357198848
	.align 4
.LC55:
	.long	2365587456
	.section	.rodata.cst8
	.align 8
.LC59:
	.long	4294967295
	.long	2146435071
	.section	.rodata.cst4
	.align 4
.LC62:
	.long	1333788672
	.align 4
.LC63:
	.long	218103808
	.align 4
.LC86:
	.long	1112014848
	.align 4
.LC87:
	.long	1096810496
	.ident	"GCC: (GNU) 4.4.7 20120313 (Red Hat 4.4.7-3)"
	.section	.note.GNU-stack,"",@progbits
