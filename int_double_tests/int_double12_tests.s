	.file	"int_double12_tests.cpp"
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC1:
	.string	"[%20.18e,%20.18e]"
	.text
.globl _Z16print_int_doubleRK10int_double
	.type	_Z16print_int_doubleRK10int_double, @function
_Z16print_int_doubleRK10int_double:
.LFB1064:
	subq	$8, %rsp
.LCFI0:
	movsd	8(%rdi), %xmm1
	movsd	.LC0(%rip), %xmm0
	xorpd	%xmm0, %xmm1
	movsd	(%rdi), %xmm0
	movl	$.LC1, %edi
	movl	$2, %eax
	call	printf
	addq	$8, %rsp
	ret
.LFE1064:
	.size	_Z16print_int_doubleRK10int_double, .-_Z16print_int_doubleRK10int_double
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC2:
	.string	"Error constructing int_double, right %20.18e < left %20.18e . Exiting.\n"
	.section	.text._ZN10int_doubleC1Edd,"axG",@progbits,_ZN10int_doubleC1Edd,comdat
	.align 2
	.weak	_ZN10int_doubleC1Edd
	.type	_ZN10int_doubleC1Edd, @function
_ZN10int_doubleC1Edd:
.LFB1063:
	subq	$8, %rsp
.LCFI1:
	movapd	%xmm1, %xmm2
	ucomisd	%xmm1, %xmm0
	jbe	.L8
	movapd	%xmm0, %xmm1
	movapd	%xmm2, %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L8:
	movsd	%xmm0, (%rdi)
	movsd	.LC0(%rip), %xmm0
	xorpd	%xmm2, %xmm0
	movsd	%xmm0, 8(%rdi)
	addq	$8, %rsp
	ret
.LFE1063:
	.size	_ZN10int_doubleC1Edd, .-_ZN10int_doubleC1Edd
	.text
	.type	_GLOBAL__I__Z16print_int_doubleRK10int_double, @function
_GLOBAL__I__Z16print_int_doubleRK10int_double:
.LFB1170:
	pushq	%rbx
.LCFI2:
	subq	$16, %rsp
.LCFI3:
	movl	$_ZStL8__ioinit, %edi
	call	_ZNSt8ios_base4InitC1Ev
	movl	$__dso_handle, %edx
	movl	$_ZStL8__ioinit, %esi
	movl	$_ZNSt8ios_base4InitD1Ev, %edi
	call	__cxa_atexit
	movabsq	$4602678819172646912, %rbx
	movq	%rbx, 8(%rsp)
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
	movl	$d_one, %eax
	movl	$0, (%rax)
	movl	$1072693248, 4(%rax)
	movl	$d_one+8, %esi
	movl	$0, (%rsi)
	movl	$-1074790400, 4(%rsi)
	movabsq	$-9223372036854775808, %rcx
	movq	%rcx, c_zero+8(%rip)
	movl	$0, %edx
	movq	%rdx, c_zero(%rip)
	movq	%rcx, c_zero+24(%rip)
	movq	%rdx, c_zero+16(%rip)
	movq	(%rax), %rax
	movq	%rax, c_one(%rip)
	movq	(%rsi), %rax
	movq	%rax, c_one+8(%rip)
	movq	%rcx, c_one+24(%rip)
	movq	%rdx, c_one+16(%rip)
	movl	$0, c_half+8(%rip)
	movl	$-1075838976, c_half+12(%rip)
	movq	%rbx, c_half(%rip)
	movq	%rcx, c_half+24(%rip)
	movq	%rdx, c_half+16(%rip)
	movsd	.LC9(%rip), %xmm1
	movsd	%xmm1, delta_int_double(%rip)
	movsd	.LC10(%rip), %xmm0
	movsd	%xmm0, delta_int_double+8(%rip)
	movsd	%xmm0, delta_int_double_neg(%rip)
	movsd	%xmm1, delta_int_double_neg+8(%rip)
	movl	$delta_blow, %edi
	call	_ZN10int_doubleC1Edd
	addq	$16, %rsp
	popq	%rbx
	ret
.LFE1170:
	.size	_GLOBAL__I__Z16print_int_doubleRK10int_double, .-_GLOBAL__I__Z16print_int_doubleRK10int_double
	.section	.ctors,"aw",@progbits
	.align 8
	.quad	_GLOBAL__I__Z16print_int_doubleRK10int_double
	.text
.globl _Z9_fpu_rnddv
	.type	_Z9_fpu_rnddv, @function
_Z9_fpu_rnddv:
.LFB1167:
	subq	$24, %rsp
.LCFI4:
	call	crlibm_init
#APP
# 1477 "../includes/int_double12.0.h" 1
	stmxcsr old__SSE_cw(%rip)
# 0 "" 2
#NO_APP
	movl	old__SSE_cw(%rip), %eax
	andl	$8127, %eax
	orb	$32, %ah
	movl	%eax, new__SSE_cw(%rip)
#APP
# 1480 "../includes/int_double12.0.h" 1
	ldmxcsr new__SSE_cw(%rip)
# 0 "" 2
#NO_APP
	movl	$0, %edx
	movq	%rdx, d_zero_zero(%rip)
	movq	%rdx, d_zero_zero+8(%rip)
	movq	%rdx, d_zero(%rip)
	movq	_nze(%rip), %rax
	movq	%rax, d_zero+8(%rip)
	movq	%rax, d_neg_zero(%rip)
	movq	%rax, d_neg_zero+8(%rip)
	movq	%rax, d_neg_neg_zero(%rip)
	movq	%rdx, d_neg_neg_zero+8(%rip)
	movsd	.LC0(%rip), %xmm1
	movsd	d_pi+8(%rip), %xmm0
	xorpd	%xmm1, %xmm0
	call	log_ru
	movsd	%xmm0, 8(%rsp)
	movsd	d_pi(%rip), %xmm0
	call	log_rd
	movapd	%xmm0, %xmm2
	ucomisd	8(%rsp), %xmm0
	jbe	.L19
	movapd	%xmm0, %xmm1
	movsd	8(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L19:
	movsd	.LC0(%rip), %xmm1
	movsd	8(%rsp), %xmm0
	xorpd	%xmm1, %xmm0
	movsd	%xmm0, d_ln_pi+8(%rip)
	movsd	%xmm2, d_ln_pi(%rip)
	movsd	d_two_pi+8(%rip), %xmm0
	xorpd	%xmm1, %xmm0
	call	log_ru
	movsd	%xmm0, 16(%rsp)
	movsd	d_two_pi(%rip), %xmm0
	call	log_rd
	movapd	%xmm0, %xmm1
	ucomisd	16(%rsp), %xmm0
	jbe	.L20
	movsd	16(%rsp), %xmm0
	movl	$.LC2, %edi
	movl	$2, %eax
	call	printf
	movl	$1, %edi
	call	exit
.L20:
	movsd	.LC0(%rip), %xmm0
	movsd	16(%rsp), %xmm2
	xorpd	%xmm2, %xmm0
	movsd	%xmm0, d_ln_two_pi+8(%rip)
	movsd	%xmm1, d_ln_two_pi(%rip)
	addq	$24, %rsp
	ret
.LFE1167:
	.size	_Z9_fpu_rnddv, .-_Z9_fpu_rnddv
	.section	.rodata.str1.1
.LC13:
	.string	"Testing...."
	.text
.globl main
	.type	main, @function
main:
.LFB1168:
	pushq	%rbx
.LCFI5:
	subq	$48, %rsp
.LCFI6:
	call	_Z9_fpu_rnddv
	movl	$0, 40(%rsp)
	movl	$-1074790400, 44(%rsp)
	movl	$0, 32(%rsp)
	movl	$1072693248, 36(%rsp)
	movl	$-18103648, 24(%rsp)
	movl	$-1179367349, 28(%rsp)
	movl	$-18103648, 16(%rsp)
	movl	$968116299, 20(%rsp)
#APP
# 384 "../includes/int_double12.0.h" 1
	movapd 32(%rsp),%XMM0
	addpd 16(%rsp),%XMM0
	movapd %XMM0,(%rsp)
	
# 0 "" 2
#NO_APP
	movq	8(%rsp), %rax
	movq	%rax, 40(%rsp)
	movq	(%rsp), %rax
	movq	%rax, 32(%rsp)
	leaq	32(%rsp), %rbx
	movq	%rbx, %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
#APP
# 446 "../includes/int_double12.0.h" 1
	movapd 16(%rsp),%xmm0
	shufpd $1,%xmm0,%xmm0
	addpd 32(%rsp),%xmm0
	movapd %xmm0,(%rsp)
	
# 0 "" 2
#NO_APP
	movq	8(%rsp), %rax
	movq	%rax, 40(%rsp)
	movq	(%rsp), %rax
	movq	%rax, 32(%rsp)
	movq	%rbx, %rdi
	call	_Z16print_int_doubleRK10int_double
	movl	$10, %edi
	call	putchar
	movl	$.LC13, %edi
	call	puts
	movl	$0, %eax
	addq	$48, %rsp
	popq	%rbx
	ret
.LFE1168:
	.size	main, .-main
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
	.local	_ZStL8__ioinit
	.comm	_ZStL8__ioinit,1,1
	.local	_ZL6d_half
	.comm	_ZL6d_half,16,16
	.weakref	_ZL20__gthrw_pthread_oncePiPFvvE,pthread_once
	.weakref	_ZL27__gthrw_pthread_getspecificj,pthread_getspecific
	.weakref	_ZL27__gthrw_pthread_setspecificjPKv,pthread_setspecific
	.weakref	_ZL22__gthrw_pthread_createPmPK14pthread_attr_tPFPvS3_ES3_,pthread_create
	.weakref	_ZL22__gthrw_pthread_cancelm,pthread_cancel
	.weakref	_ZL26__gthrw_pthread_mutex_lockP15pthread_mutex_t,pthread_mutex_lock
	.weakref	_ZL29__gthrw_pthread_mutex_trylockP15pthread_mutex_t,pthread_mutex_trylock
	.weakref	_ZL28__gthrw_pthread_mutex_unlockP15pthread_mutex_t,pthread_mutex_unlock
	.weakref	_ZL26__gthrw_pthread_mutex_initP15pthread_mutex_tPK19pthread_mutexattr_t,pthread_mutex_init
	.weakref	_ZL30__gthrw_pthread_cond_broadcastP14pthread_cond_t,pthread_cond_broadcast
	.weakref	_ZL25__gthrw_pthread_cond_waitP14pthread_cond_tP15pthread_mutex_t,pthread_cond_wait
	.weakref	_ZL26__gthrw_pthread_key_createPjPFvPvE,pthread_key_create
	.weakref	_ZL26__gthrw_pthread_key_deletej,pthread_key_delete
	.weakref	_ZL30__gthrw_pthread_mutexattr_initP19pthread_mutexattr_t,pthread_mutexattr_init
	.weakref	_ZL33__gthrw_pthread_mutexattr_settypeP19pthread_mutexattr_ti,pthread_mutexattr_settype
	.weakref	_ZL33__gthrw_pthread_mutexattr_destroyP19pthread_mutexattr_t,pthread_mutexattr_destroy
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC0:
	.long	0
	.long	-2147483648
	.long	0
	.long	0
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC9:
	.long	0
	.long	1048576
	.align 8
.LC10:
	.long	0
	.long	-2146435072
	.section	.eh_frame,"a",@progbits
.Lframe1:
	.long	.LECIE1-.LSCIE1
.LSCIE1:
	.long	0x0
	.byte	0x1
.globl __gxx_personality_v0
	.string	"zPR"
	.uleb128 0x1
	.sleb128 -8
	.byte	0x10
	.uleb128 0x6
	.byte	0x3
	.long	__gxx_personality_v0
	.byte	0x3
	.byte	0xc
	.uleb128 0x7
	.uleb128 0x8
	.byte	0x90
	.uleb128 0x1
	.align 8
.LECIE1:
.LSFDE1:
	.long	.LEFDE1-.LASFDE1
.LASFDE1:
	.long	.LASFDE1-.Lframe1
	.long	.LFB1064
	.long	.LFE1064-.LFB1064
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI0-.LFB1064
	.byte	0xe
	.uleb128 0x10
	.align 8
.LEFDE1:
.LSFDE3:
	.long	.LEFDE3-.LASFDE3
.LASFDE3:
	.long	.LASFDE3-.Lframe1
	.long	.LFB1063
	.long	.LFE1063-.LFB1063
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI1-.LFB1063
	.byte	0xe
	.uleb128 0x10
	.align 8
.LEFDE3:
.LSFDE5:
	.long	.LEFDE5-.LASFDE5
.LASFDE5:
	.long	.LASFDE5-.Lframe1
	.long	.LFB1170
	.long	.LFE1170-.LFB1170
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI2-.LFB1170
	.byte	0xe
	.uleb128 0x10
	.byte	0x4
	.long	.LCFI3-.LCFI2
	.byte	0xe
	.uleb128 0x20
	.byte	0x83
	.uleb128 0x2
	.align 8
.LEFDE5:
.LSFDE7:
	.long	.LEFDE7-.LASFDE7
.LASFDE7:
	.long	.LASFDE7-.Lframe1
	.long	.LFB1167
	.long	.LFE1167-.LFB1167
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI4-.LFB1167
	.byte	0xe
	.uleb128 0x20
	.align 8
.LEFDE7:
.LSFDE9:
	.long	.LEFDE9-.LASFDE9
.LASFDE9:
	.long	.LASFDE9-.Lframe1
	.long	.LFB1168
	.long	.LFE1168-.LFB1168
	.uleb128 0x0
	.byte	0x4
	.long	.LCFI5-.LFB1168
	.byte	0xe
	.uleb128 0x10
	.byte	0x4
	.long	.LCFI6-.LCFI5
	.byte	0xe
	.uleb128 0x40
	.byte	0x83
	.uleb128 0x2
	.align 8
.LEFDE9:
	.ident	"GCC: (GNU) 4.3.3"
	.section	.note.GNU-stack,"",@progbits
